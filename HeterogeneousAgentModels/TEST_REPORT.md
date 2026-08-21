# Heterogeneous-agent model test report

Retested on 21 August 2026 with MATLAB R2026a Update 4 and one NVIDIA GeForce RTX 5080 Laptop GPU. The examples repository was at commit `59fd543f5b91f2c9f0b120b60f966dd04db0152a`, and VFI Toolkit was at commit `b97fe98ea93101b9603e57fd43ee51cde0eb0900`.

The comparison results are from the 13 and 19 August tests. Those used examples commits `fbecce303f6a95b7a03aa3ace918cd2c29031877` and `38341f855639690143087714d2dc584fcbb789d1`, and VFI Toolkit commits `eb44124b0d68a5d6f2155480fdf4377aec1fa64b` and `b0ddb078def8be0b9b6cdd69db9033bc38f5df49`.

Each run used a fresh MATLAB process and the separately checked-out `VFIToolkit-matlab` repository on the MATLAB path. Tests ran from temporary directories so `GL2017A.mat` and `GL2017B.mat` were not written to the examples repository. Both modes used `vfoptions.ngridinterp=15`, with matching `simoptions` fields. No model grids, parameters, tolerances, or other model options were changed. Timing instrumentation was inserted only into the in-memory script text used by the temporary test driver.

## Results

| Script | `gridinterplayer` | Previous runtime | 21 August runtime | Runtime change | Convergence evidence and result |
|---|---:|---:|---:|---:|---|
| `Aiyagari1994.m` | 0 | 15.4 s | 12.3 s | 20% faster | Pass; finite stationary GE and capital-market residual `-8.52e-5` |
| `Aiyagari1994.m` | 1 | 26.4 s | 6.3 s | 76% faster | Pass; finite stationary GE and capital-market residual `-1.32e-4` |
| `Aiyagari1994TransitionPath.m` | 0 | 89.7 s | 81.4 s | 9% faster | Pass; converged in 29 iterations with final ratio `0.87` |
| `Aiyagari1994TransitionPath.m` | 1 | 195.9 s | 93.5 s | 52% faster | Pass; converged in 29 iterations with final ratio `0.88` |
| `GuerrieriLorenzoni2017_Example.m` | 0 | 2,274.3 s | 1,913.5 s | 16% faster | Pass; flexible and sticky transitions converged in 310 and 30 iterations, both with displayed final ratios `1.00` |
| `GuerrieriLorenzoni2017_Example.m` | 1 | 3,093 s before stop at sticky iteration 121 | 2,753.1 s before requested stop at the start of sticky iteration 90 | Not directly comparable | Partial pass; flexible transition and Figure 6 completed, but sticky transition again diverged; last completed ratio `16.66` at iteration 89 |

All five directly comparable completed runs were faster with VFI Toolkit commit `b97fe98e`. The largest gains were in grid-interpolation mode for the stationary and transition Aiyagari examples. The full mode-0 Guerrieri--Lorenzoni run improved from 37.9 to 31.9 minutes, but it remains well above the requested ten-minute threshold.

The mode-1 Guerrieri--Lorenzoni flexible transition converged in 309 iterations with a final displayed ratio of `1.00`. The transition call itself took 1,656.5 seconds (27.6 minutes). The nearest previous measure was a 31.2-minute checkpoint after the flexible transition and subsequent value, policy, and distribution calculations, so it suggests a modest improvement but is not an exact like-for-like timing. Figure 6 completed successfully again. The sticky transition began 2,259.1 seconds (37.7 minutes) into the run.

## Sticky-transition result

The mode-1 sticky transition still does **not** converge under the unchanged example settings. Its ratio initially fell from `17.50` at iteration 2 to a displayed `1.00` at iterations 15--23 without crossing the stopping threshold. It then turned upward: `1.90` at iteration 43, `5.65` at iteration 67, `8.26` at iteration 75, and `16.66` at iteration 89. These checkpoints reproduce the previous divergent path exactly through iteration 89.

The user requested termination while the divergence was being monitored. The MATLAB output stream closed as the process was stopped at the start of iteration 90; that output-stream error was a termination artifact, not a model or toolkit failure. Continuing toward the default 1,000-iteration cap was not warranted because the convergence distance was increasing rapidly.

For clarity, the mode-0 sticky transition does converge: it completed in 30 iterations and its transition call took 119.7 seconds. The failure is specific to the grid-interpolation (`gridinterplayer=1`) sticky transition under the example's current settings.

## Grid-interpolation status

The earlier `SimPanelValues_TransPath_InfHorz` indexing failure remains fixed. All four Figure 6 simulation calls completed in the latest mode-1 run. The unresolved problem is the sticky-transition convergence path, not the former panel-simulation indexing error.
