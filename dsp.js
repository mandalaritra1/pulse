// Pulse — signal processing core.
// Camera PPG: the fingertip's brightness dips slightly with every heartbeat.
// We resample the per-frame brightness onto a uniform grid, strip the slow
// drift, smooth out sensor noise, and find the dominant period by
// autocorrelation (robust), cross-checked with explicit peak detection.
(function (root) {
  'use strict';

  const FS = 30; // uniform resample rate, Hz

  // Linear-interpolate irregular (t [s], v) samples onto a uniform grid
  // covering the last `dur` seconds ending at the newest sample.
  function resample(ts, vs, dur, fs = FS) {
    const n = ts.length;
    if (n < 2) return [];
    const tEnd = ts[n - 1];
    const tStart = Math.max(ts[0], tEnd - dur);
    const m = Math.floor((tEnd - tStart) * fs) + 1;
    const out = new Array(m);
    let j = 0;
    while (j < n - 1 && ts[j + 1] < tStart) j++;
    for (let i = 0; i < m; i++) {
      const t = tStart + i / fs;
      while (j < n - 2 && ts[j + 1] < t) j++;
      const t0 = ts[j], t1 = ts[j + 1];
      const a = t1 > t0 ? Math.min(1, Math.max(0, (t - t0) / (t1 - t0))) : 0;
      out[i] = vs[j] + a * (vs[j + 1] - vs[j]);
    }
    return out;
  }

  // Centred moving average with half-width h (window shrinks at the edges).
  function movingAverage(x, h) {
    const n = x.length, out = new Array(n);
    const cs = new Float64Array(n + 1);
    for (let i = 0; i < n; i++) cs[i + 1] = cs[i] + x[i];
    for (let i = 0; i < n; i++) {
      const a = Math.max(0, i - h), b = Math.min(n, i + h + 1);
      out[i] = (cs[b] - cs[a]) / (b - a);
    }
    return out;
  }

  // Band-pass by difference of moving averages:
  //   high-pass: subtract ~1.5 s trend (kills drift / exposure changes,
  //              still passes 40 bpm)
  //   low-pass:  ~0.17 s smoothing     (kills frame noise)
  // Sign is flipped so systolic beats point UP.
  function bandpass(x, fs = FS, smoothHalf = Math.max(1, Math.round(0.05 * fs))) {
    if (x.length < 5) return x.slice();
    const trend = movingAverage(x, Math.round(0.75 * fs));
    const hp = x.map((v, i) => trend[i] - v);
    return smoothHalf > 0 ? movingAverage(hp, smoothHalf) : hp;
  }

  function standardize(x) {
    const n = x.length;
    let m = 0; for (const v of x) m += v; m /= n;
    let s = 0; for (const v of x) s += (v - m) * (v - m);
    s = Math.sqrt(s / n) || 1;
    return x.map(v => (v - m) / s);
  }

  // Dominant beat period by normalised autocorrelation.
  // Returns { bpm, q } where q in [0,1] is the correlation at that lag
  // (a clean pulse gives q > 0.6; noise gives q < 0.3).
  function autocorrBpm(sig, fs = FS, minBpm = 35, maxBpm = 210) {
    const x = standardize(sig);
    const n = x.length;
    const lagMin = Math.floor(fs * 60 / maxBpm);
    const lagMax = Math.min(Math.ceil(fs * 60 / minBpm), Math.floor(n * 0.6));
    if (lagMax <= lagMin + 2) return { bpm: 0, q: 0 };
    const Lall = Math.floor(n * 0.6);
    const r = new Float64Array(Lall + 2);
    for (let L = 1; L <= Lall + 1 && L < n; L++) {
      let s = 0;
      for (let i = 0; i + L < n; i++) s += x[i] * x[i + L];
      r[L] = s / (n - L);
    }
    // local maxima in range
    const peaks = [];
    for (let L = lagMin; L <= lagMax; L++) {
      if (r[L] > r[L - 1] && r[L] >= r[L + 1] && r[L] > 0) peaks.push(L);
    }
    if (!peaks.length) return { bpm: 0, q: 0 };
    // Subharmonic summation: the true period L also correlates at 2L, 3L...
    // while a multiple (2T, 3T) or half-period (T/2) does not collect that
    // support. Score each candidate by the sum over its multiples.
    const rAt = (t) => {
      const i = Math.floor(t), f = t - i;
      return i + 1 <= Lall ? r[i] * (1 - f) + r[i + 1] * f : 0;
    };
    let best = peaks[0], bestScore = -Infinity;
    for (const L of peaks) {
      let s = 0;
      for (let k = 1; k * L <= Lall; k++) {
        // allow a little smear at long lags (heart-rate variability)
        const c = k * L, w = Math.max(1, Math.round(0.04 * c));
        let m = -Infinity;
        for (let d = -w; d <= w; d++) m = Math.max(m, rAt(c + d));
        s += m;
      }
      if (s > bestScore) { bestScore = s; best = L; }
    }
    // parabolic interpolation for sub-sample lag
    const y0 = r[best - 1], y1 = r[best], y2 = r[best + 1];
    const den = y0 - 2 * y1 + y2;
    const d = den !== 0 ? Math.max(-0.5, Math.min(0.5, 0.5 * (y0 - y2) / den)) : 0;
    const lag = best + d;
    return { bpm: 60 * fs / lag, q: Math.max(0, Math.min(1, y1)) };
  }

  // Beat peaks: local maxima above an adaptive threshold, with a refractory
  // period tied to the expected beat interval when known.
  function findPeaks(sig, fs = FS, expectedBpm = 0) {
    const n = sig.length;
    if (n < 5) return [];
    // robust threshold: a fraction of the 90th percentile, so one noise
    // spike can't hide the real beats
    const sorted = sig.slice().sort((a, b) => a - b);
    const thr = 0.3 * sorted[Math.floor(0.9 * (n - 1))];
    const refr = expectedBpm > 0 ? 0.6 * 60 / expectedBpm : 60 / 200;
    const minGap = Math.round(refr * fs);
    const out = [];
    for (let i = 2; i < n - 2; i++) {
      const v = sig[i];
      if (v > thr && v >= sig[i - 1] && v >= sig[i - 2] && v > sig[i + 1] && v >= sig[i + 2]) {
        if (out.length && i - out[out.length - 1] < minGap) {
          if (v > sig[out[out.length - 1]]) out[out.length - 1] = i;
        } else out.push(i);
      }
    }
    return out;
  }

  // BPM from median inter-beat interval.
  function peakBpm(peaks, fs = FS) {
    if (peaks.length < 3) return 0;
    const ibi = [];
    for (let i = 1; i < peaks.length; i++) ibi.push((peaks[i] - peaks[i - 1]) / fs);
    ibi.sort((a, b) => a - b);
    const med = ibi[Math.floor(ibi.length / 2)];
    return 60 / med;
  }

  // Full estimate on the last `dur` seconds of raw (t, value) samples.
  function estimate(ts, vs, dur = 8, fs = FS) {
    const u = resample(ts, vs, dur, fs);
    if (u.length < fs * 3) return { bpm: 0, q: 0, sig: [], peaks: [] };
    const sig = bandpass(u, fs);
    const ac = autocorrBpm(sig, fs);
    const peaks = findPeaks(sig, fs, ac.bpm);
    const pb = peakBpm(peaks, fs);
    // agreement between the two methods boosts/penalises confidence
    let q = ac.q;
    if (pb > 0 && ac.bpm > 0) {
      const rel = Math.abs(pb - ac.bpm) / ac.bpm;
      if (rel > 0.25) q *= 0.75;
    }
    return { bpm: ac.bpm, q, peakBpm: pb, sig, peaks };
  }

  function median(a) {
    if (!a.length) return 0;
    const s = a.slice().sort((x, y) => x - y);
    const m = s.length >> 1;
    return s.length % 2 ? s[m] : 0.5 * (s[m - 1] + s[m]);
  }

  // ======================================================================
  //  Beat-to-beat timing (for 30 s and HRV modes)
  // ======================================================================
  // HRV needs beat times to a few ms, but frames are 33 ms apart. We time
  // each beat at the point of steepest upstroke (the most stable fiducial
  // on a PPG pulse) and refine it to sub-frame precision by fitting a
  // parabola through the slope around its maximum.

  // Resample samples [i0, i1) onto a uniform grid. Returns {t0, x}.
  function resampleRange(ts, vs, i0, i1, fs = FS) {
    const tA = ts[i0], tB = ts[i1 - 1];
    const m = Math.floor((tB - tA) * fs) + 1;
    const x = new Array(m);
    let j = i0;
    for (let i = 0; i < m; i++) {
      const t = tA + i / fs;
      while (j < i1 - 2 && ts[j + 1] < t) j++;
      const dt = ts[j + 1] - ts[j];
      const a = dt > 0 ? Math.min(1, Math.max(0, (t - ts[j]) / dt)) : 0;
      x[i] = vs[j] + a * (vs[j + 1] - vs[j]);
    }
    return { t0: tA, x };
  }

  // Beats in all segments. `sg` labels each sample with a segment id; a new
  // segment starts whenever the finger was lifted, and no interval is ever
  // measured across a gap. Returns [{t, seg}] in seconds.
  function detectBeats(ts, vs, sg, fs = FS) {
    const beats = [];
    let i0 = 0;
    while (i0 < ts.length) {
      let i1 = i0;
      while (i1 < ts.length && sg[i1] === sg[i0]) i1++;
      if (i1 - i0 > 10 && ts[i1 - 1] - ts[i0] >= 5) beatsInSegment(ts, vs, i0, i1, sg[i0], fs, beats);
      i0 = i1;
    }
    return beats;
  }

  function beatsInSegment(ts, vs, i0, i1, seg, fs, out) {
    const { t0, x } = resampleRange(ts, vs, i0, i1, fs);
    const s = bandpass(x, fs, 1); // light smoothing keeps the upstroke sharp
    const n = s.length;
    // expected beat period, re-estimated every 10 s (HR drifts in long runs)
    const B = 10 * fs, per = new Float64Array(n);
    for (let b = 0; b < n; b += B) {
      const a = Math.max(0, Math.min(b, n - B)), e = Math.min(n, a + B);
      const ac = autocorrBpm(s.slice(a, e), fs);
      const p = ac.bpm > 0 ? 60 / ac.bpm : 0.8;
      for (let i = b; i < Math.min(n, b + B); i++) per[i] = p;
    }
    // local envelope (moving max over +-1 s) for an adaptive threshold
    const h = fs, env = new Float64Array(n);
    for (let i = 0; i < n; i++) {
      let m = 0;
      for (let k = Math.max(0, i - h); k <= Math.min(n - 1, i + h); k++) m = Math.max(m, s[k]);
      env[i] = m;
    }
    // systolic peaks with a refractory period of ~half a beat
    const pk = [];
    for (let i = 2; i < n - 2; i++) {
      const v = s[i];
      if (v > 0.35 * env[i] && v >= s[i - 1] && v >= s[i - 2] && v > s[i + 1] && v >= s[i + 2]) {
        const gap = 0.55 * per[i] * fs;
        if (pk.length && i - pk[pk.length - 1] < gap) {
          if (v > s[pk[pk.length - 1]]) pk[pk.length - 1] = i;
        } else pk.push(i);
      }
    }
    // coarse fiducial: maximum slope on the upstroke before each peak
    const d = (k) => (k > 0 && k < n - 1 ? s[k + 1] - s[k - 1] : 0);
    let bt = [];
    for (const p of pk) {
      const lo = Math.max(1, p - Math.round(0.45 * per[p] * fs));
      let k = p, best = -Infinity;
      for (let j = lo; j <= p; j++) if (d(j) > best) { best = d(j); k = j; }
      if (best <= 0) continue;
      const y0 = d(k - 1), y1 = d(k), y2 = d(k + 1), den = y0 - 2 * y1 + y2;
      const delta = den !== 0 ? Math.max(-0.5, Math.min(0.5, 0.5 * (y0 - y2) / den)) : 0;
      bt.push(t0 + (k + delta) / fs);
    }
    if (bt.length < 4) { for (const t of bt) out.push({ t, seg, q: 0.5 }); return; }

    // ---- fine timing by template matching on the RAW samples ----
    // High-passed raw values (no interpolation of the data itself).
    const trend = movingAverage(x, Math.round(0.75 * fs));
    const trendAt = (t) => {
      const u = (t - t0) * fs, i = Math.max(0, Math.min(n - 2, Math.floor(u))), a = Math.max(0, Math.min(1, u - i));
      return trend[i] + a * (trend[i + 1] - trend[i]);
    };
    const tsS = ts.slice(i0, i1), y = new Float64Array(i1 - i0);
    for (let j = 0; j < y.length; j++) y[j] = trendAt(tsS[j]) - vs[i0 + j];
    // window around each beat, kept inside one beat at high heart rates
    const ibisC = []; for (let i = 1; i < bt.length; i++) ibisC.push(bt[i] - bt[i - 1]);
    const P = median(ibisC);
    const WA = Math.min(0.15, 0.35 * P), WB = Math.max(0.08, Math.min(0.3, 0.75 * P - WA));
    const BIN = 0.005, NB = Math.round((WA + WB) / BIN) + 1;
    // first sample index with ts >= t (binary search)
    const lowerBound = (t) => { let a = 0, b = tsS.length; while (a < b) { const m = (a + b) >> 1; if (tsS[m] < t) a = m + 1; else b = m; } return a; };

    let q = new Array(bt.length).fill(1);
    const sig = new Array(bt.length).fill(0);
    for (let iter = 0; iter < 2; iter++) {
      // Template: frames fall at random phase relative to each beat, so
      // stacking all beats samples the pulse shape every 5 ms
      // ("equivalent-time sampling").
      const sum = new Float64Array(NB), cnt = new Float64Array(NB);
      bt.forEach((tb, k) => {
        if (q[k] < 0.5) return;
        for (let j = lowerBound(tb - WA); j < tsS.length && tsS[j] <= tb + WB; j++) {
          const b = Math.round((tsS[j] - tb + WA) / BIN);
          if (b >= 0 && b < NB) { sum[b] += y[j]; cnt[b]++; }
        }
      });
      let T = Array.from(sum, (v, i) => (cnt[i] ? v / cnt[i] : 0));
      T = movingAverage(T, 2);
      const Tm = mean(T);
      T = T.map((v, i) => (v - Tm) * (0.5 - 0.5 * Math.cos(2 * Math.PI * i / (NB - 1)))); // taper edges
      const tmpl = (tau) => {
        const u = (tau + WA) / BIN, i = Math.floor(u);
        if (i < 0 || i >= NB - 1) return 0;
        return T[i] + (u - i) * (T[i + 1] - T[i]);
      };
      // Fit each beat as  y_j = a*T(t_j - tb - d) + c  by least squares over
      // a FIXED set of the samples that actually exist (so dropped frames
      // just lose information instead of biasing the timing). For each
      // shift d the best a, c are analytic; we maximise the correlation.
      const STEP = 0.002, R = 35, SH = R * STEP;
      bt = bt.map((tb, k) => {
        const jA = lowerBound(tb - WA - SH);
        let jB = jA; while (jB < tsS.length && tsS[jB] <= tb + WB + SH) jB++;
        const nS = jB - jA;
        if (nS < 6) { q[k] = 0; return tb; }
        let Sy = 0, Syy = 0;
        for (let j = jA; j < jB; j++) { Sy += y[j]; Syy += y[j] * y[j]; }
        const vy = Syy - Sy * Sy / nS;
        const fit = (dd) => { // correlation of samples with the shifted template
          let St = 0, Stt = 0, Sty = 0;
          for (let j = jA; j < jB; j++) { const T0 = tmpl(tsS[j] - tb - dd); St += T0; Stt += T0 * T0; Sty += T0 * y[j]; }
          const cov = Sty - St * Sy / nS, vt = Stt - St * St / nS;
          return vt > 0 && vy > 0 ? cov / Math.sqrt(vt * vy) : -1;
        };
        let bestD = 0, bestS = -Infinity; const sc = [];
        for (let m = -R; m <= R; m++) { const v = fit(m * STEP); sc.push(v); if (v > bestS) { bestS = v; bestD = m; } }
        const c = bestD + R;
        let frac = 0;
        if (c > 0 && c < sc.length - 1) {
          const den = sc[c - 1] - 2 * sc[c] + sc[c + 1];
          frac = den !== 0 ? Math.max(-0.5, Math.min(0.5, 0.5 * (sc[c - 1] - sc[c + 1]) / den)) : 0;
        }
        q[k] = bestS; // per-beat quality: correlation with the average pulse
        const tn = tb + (bestD + frac) * STEP;
        // timing uncertainty of this beat (Cramer-Rao style):
        //   sigma_t^2 = sigma_noise^2 / (a^2 * sum (T'(t_j) - mean T')^2)
        if (iter === 1) {
          const dsh = tn - tb, h = 0.004;
          let St = 0, Stt = 0, Sty = 0, Sd = 0, Sdd = 0;
          for (let j = jA; j < jB; j++) {
            const tau = tsS[j] - tb - dsh, T0 = tmpl(tau);
            const dT = (tmpl(tau + h) - tmpl(tau - h)) / (2 * h);
            St += T0; Stt += T0 * T0; Sty += T0 * y[j]; Sd += dT; Sdd += dT * dT;
          }
          const vt = Stt - St * St / nS, cov = Sty - St * Sy / nS;
          const a = vt > 0 ? cov / vt : 0;
          const resVar = Math.max(0, (vy - (vt > 0 ? cov * cov / vt : 0)) / Math.max(1, nS - 3));
          const info = a * a * (Sdd - Sd * Sd / nS);
          sig[k] = info > 0 ? Math.sqrt(resVar / info) : 0.05;
        }
        return tn;
      });
    }
    bt.forEach((t, k) => out.push({ t, seg, q: q[k], sig: sig[k] }));
  }

  // Inter-beat intervals (ms) with artifact flags: physiologically
  // impossible, or >20% away from the local median (missed / extra beat).
  function intervals(beats) {
    const out = [];
    for (let i = 1; i < beats.length; i++) {
      if (beats[i].seg === beats[i - 1].seg) {
        const shapeOk = (beats[i].q ?? 1) >= 0.5 && (beats[i - 1].q ?? 1) >= 0.5;
        out.push({ t: beats[i].t, ibi: 1000 * (beats[i].t - beats[i - 1].t), seg: beats[i].seg, valid: shapeOk,
          sig: 1000 * (beats[i].sig || 0) }); // timing sd of the closing beat, ms
      }
    }
    for (let i = 0; i < out.length; i++) {
      const v = out[i].ibi, win = [];
      for (let k = i - 5; k <= i + 5; k++) if (k >= 0 && k < out.length && k !== i) win.push(out[k].ibi);
      const med = median(win.length ? win : [v]);
      out[i].valid = out[i].valid && v >= 270 && v <= 2000 && Math.abs(v - med) <= 0.2 * med;
    }
    return out;
  }

  function mean(a) { let s = 0; for (const v of a) s += v; return a.length ? s / a.length : 0; }
  function sd(a) {
    if (a.length < 2) return 0;
    const m = mean(a); let s = 0;
    for (const v of a) s += (v - m) * (v - m);
    return Math.sqrt(s / (a.length - 1));
  }
  function pct(sorted, p) { return sorted[Math.min(sorted.length - 1, Math.max(0, Math.round(p * (sorted.length - 1))))]; }

  // Lomb-Scargle periodogram: handles the uneven beat-to-beat sampling and
  // the gaps from lifted fingers directly, no interpolation needed.
  function lombScargle(t, y, freqs) {
    return freqs.map((f) => {
      const w = 2 * Math.PI * f;
      let s2 = 0, c2 = 0;
      for (const ti of t) { s2 += Math.sin(2 * w * ti); c2 += Math.cos(2 * w * ti); }
      const tau = Math.atan2(s2, c2) / (2 * w);
      let yc = 0, ys = 0, cc = 0, ss = 0;
      for (let i = 0; i < t.length; i++) {
        const c = Math.cos(w * (t[i] - tau)), s = Math.sin(w * (t[i] - tau));
        yc += y[i] * c; ys += y[i] * s; cc += c * c; ss += s * s;
      }
      return 0.5 * ((cc > 0 ? yc * yc / cc : 0) + (ss > 0 ? ys * ys / ss : 0));
    });
  }

  // Everything HRV. `ibis` from intervals().
  function hrvStats(ibis) {
    const nnI = ibis.filter((x) => x.valid);
    const nn = nnI.map((x) => x.ibi);
    const N = nn.length, total = ibis.length;
    const artifactPct = total ? 100 * (1 - N / total) : 100;
    if (N < 10) return { ok: false, n: N, total, artifactPct };

    const meanNN = mean(nn), sdnn = sd(nn);
    // successive differences only between adjacent, clean intervals
    const dd = [];
    for (let i = 0; i + 1 < ibis.length; i++) {
      const a = ibis[i], b = ibis[i + 1];
      if (a.valid && b.valid && a.seg === b.seg) dd.push(b.ibi - a.ibi);
    }
    const rmssdRaw = Math.sqrt(mean(dd.map((v) => v * v)));
    const pnn50 = dd.length ? 100 * dd.filter((v) => Math.abs(v) > 50).length / dd.length : 0;
    // Beat-timing noise (sd s per beat) adds 6 s^2 to RMSSD^2 and 2 s^2 to
    // SDNN^2. Subtract that noise floor, using each beat's own uncertainty.
    const s2 = mean(nnI.map((x) => x.sig * x.sig));
    const timingSd = Math.sqrt(s2);
    const rmssd = Math.sqrt(Math.max(0, rmssdRaw * rmssdRaw - 6 * s2));
    const sdnnC = Math.sqrt(Math.max(0, sdnn * sdnn - 2 * s2));
    const sd1Raw = Math.sqrt(0.5) * sd(dd);
    const sd1 = Math.sqrt(Math.max(0, sd1Raw * sd1Raw - 3 * s2));
    const sd2 = Math.sqrt(Math.max(0, 2 * sdnnC * sdnnC - sd1 * sd1));

    // heart rate: mean, and range of a 5-beat running median
    const hr = nn.map((v) => 60000 / v);
    const hrSm = hr.map((_, i) => median(hr.slice(Math.max(0, i - 2), i + 3)));
    const meanHR = 60000 / meanNN;

    // Baevsky stress index (50 ms histogram)
    const sorted = nn.slice().sort((a, b) => a - b);
    const bins = new Map();
    for (const v of nn) { const b = Math.floor(v / 50); bins.set(b, (bins.get(b) || 0) + 1); }
    let modeBin = 0, modeCount = 0;
    for (const [b, c] of bins) if (c > modeCount) { modeCount = c; modeBin = b; }
    const Mo = (modeBin * 50 + 25) / 1000, AMo = 100 * modeCount / N;
    const MxDMn = Math.max(0.05, (pct(sorted, 0.975) - pct(sorted, 0.025)) / 1000);
    const stressIndex = AMo / (2 * Mo * MxDMn);

    const span = nnI[nnI.length - 1].t - nnI[0].t;
    const res = {
      ok: true, n: N, total, artifactPct, span,
      meanHR, minHR: Math.min(...hrSm), maxHR: Math.max(...hrSm),
      meanNN, sdnn: sdnnC, rmssd, pnn50, sd1, sd2, stressIndex,
      sdnnRaw: sdnn, rmssdRaw, timingSd,
    };

    // frequency domain (needs >= ~2 min for LF to mean anything)
    if (span >= 50 && N >= 40) {
      const tt = nnI.map((x) => x.t);
      // remove linear trend
      const tm = mean(tt), ym = mean(nn);
      let sxy = 0, sxx = 0;
      for (let i = 0; i < N; i++) { sxy += (tt[i] - tm) * (nn[i] - ym); sxx += (tt[i] - tm) ** 2; }
      const slope = sxx > 0 ? sxy / sxx : 0;
      const y = nn.map((v, i) => v - ym - slope * (tt[i] - tm));
      const df = 0.002, freqs = [];
      for (let f = 0.004; f <= 0.5; f += df) freqs.push(f);
      const P = lombScargle(tt, y, freqs);
      // scale so the spectrum integrates to the (detrended) variance
      const varY = mean(y.map((v) => v * v));
      const tot = P.reduce((s, v) => s + v, 0) * df || 1;
      const psd = P.map((v) => v * varY / tot); // ms^2/Hz
      const band = (a, b) => freqs.reduce((s, f, i) => (f >= a && f < b ? s + psd[i] * df : s), 0);
      // breathing: dominant peak 0.1-0.5 Hz (6-30 breaths/min)
      let fb = 0, pb = 0;
      freqs.forEach((f, i) => { if (f >= 0.1 && f <= 0.5 && psd[i] > pb) { pb = psd[i]; fb = f; } });
      const inBand = psd.filter((_, i) => freqs[i] >= 0.1 && freqs[i] <= 0.5).sort((a, b) => a - b);
      res.respRate = pb > 3 * pct(inBand, 0.5) ? 60 * fb : null;
      res.spectrum = { freqs, psd };
      if (span >= 110) {
        const lf = band(0.04, 0.15), hf = band(0.15, 0.4);
        Object.assign(res, { lf, hf, lfhf: hf > 0 ? lf / hf : null, lfnu: 100 * lf / (lf + hf), hfnu: 100 * hf / (lf + hf) });
      }
    }
    return res;
  }

  const api = {
    FS, resample, movingAverage, bandpass, standardize, autocorrBpm, findPeaks, peakBpm, estimate, median,
    resampleRange, detectBeats, intervals, hrvStats, lombScargle, mean, sd,
  };
  if (typeof module !== 'undefined' && module.exports) module.exports = api;
  else root.PulseDSP = api;
})(typeof self !== 'undefined' ? self : this);
