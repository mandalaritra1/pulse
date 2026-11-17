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
  function bandpass(x, fs = FS) {
    if (x.length < 5) return x.slice();
    const trend = movingAverage(x, Math.round(0.75 * fs));
    const hp = x.map((v, i) => trend[i] - v);
    return movingAverage(hp, Math.max(1, Math.round(0.05 * fs)));
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

  const api = { FS, resample, movingAverage, bandpass, standardize, autocorrBpm, findPeaks, peakBpm, estimate, median };
  if (typeof module !== 'undefined' && module.exports) module.exports = api;
  else root.PulseDSP = api;
})(typeof self !== 'undefined' ? self : this);
