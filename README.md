# Pulse

Ad-free camera heart-rate monitor (PPG). Static web app, no server, no tracking.

## Put it on your phone (GitHub Pages, ~3 min)

1. On github.com create a new **public** repo, e.g. `pulse`.
2. **Add file → Upload files** and drop in everything in this folder.
3. **Settings → Pages →** Source: *Deploy from a branch*, branch `main`, folder `/ (root)` → Save.
4. After ~1 minute open `https://<your-username>.github.io/pulse/` in **Chrome on Android**.
5. Chrome menu (⋮) → **Add to home screen / Install app**. It now opens full-screen and works offline.

(The camera needs HTTPS, which is why it can't run from a local file or inside a claude.ai artifact.)

## Files

- `index.html` – UI, camera + torch control, measurement logic
- `dsp.js` – signal processing (resample → band-pass → autocorrelation with
  subharmonic summation, cross-checked by peak detection)
- `sw.js`, `manifest.webmanifest`, `icon-*.png` – offline + installable

## Modes

- **Quick**: about 8 s; stops as soon as the estimate is stable.
- **30 s**: averages over several breaths; best for resting heart rate. Also gives range and a 30 s RMSSD.
- **HRV** (1 / 2 / 5 min): RMSSD, SDNN, pNN50, SD1/SD2, Baevsky stress index, breathing rate,
  LF/HF (≥ 2 min), tachogram, Poincaré plot, spectrum, CSV export of every beat.
  Lifting the finger pauses the recording; intervals across a gap are never used.

### How HRV timing works

Beats are found on the band-passed signal, then each one is timed by least-squares fitting an
average pulse template (built by stacking all beats, which samples the pulse shape every 5 ms
because frames land at random phases) to the raw samples, with free amplitude and offset.
Dropped frames therefore lose information but don't bias timing. Each beat also gets a
Cramér–Rao-style timing uncertainty; its contribution (6σ² to RMSSD², 2σ² to SDNN²) is subtracted.
Spectra use Lomb–Scargle on the uneven beat series. Validated on synthetic data with known
ground truth (`tests/` in the source session): typically within ~5% on all metrics.

## How a quick reading finishes

Estimates every 0.25 s on the last ≤8 s of data. The result is reported as soon as
2.5 s of consecutive estimates agree within max(3 bpm, 4%), all with clean signal, after
at least 6.5 s (and 7 beats). Typically 7–9 s after your finger goes on. Gives up at 25 s.

If you change any file, bump `VERSION` in `sw.js` and the `?v=` on `dsp.js` (in index.html and sw.js) so installed copies update.
