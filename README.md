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

## How a reading finishes

Estimates every 0.25 s on the last ≤8 s of data. The result is reported as soon as
2.5 s of consecutive estimates agree within max(3 bpm, 4%), all with clean signal, after
at least 6.5 s (and 7 beats). Typically 7–9 s after your finger goes on. Gives up at 25 s.

If you change any file, bump `VERSION` in `sw.js` so installed copies update.
