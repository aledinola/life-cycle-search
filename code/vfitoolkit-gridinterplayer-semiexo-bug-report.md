# Possible VFI Toolkit Bug With `gridinterplayer = 1` and Semi-Exogenous States

I found what looks like a bug in the finite-horizon stationary-distribution code for models with semi-exogenous states when `gridinterplayer = 1`.

Running `main.m` in this repository as committed reproduces the bug.

## Summary

- The model runs through value function iteration with `vfoptions.gridinterplayer = 0`.
- The same model fails when `vfoptions.gridinterplayer = 1`.
- The failure happens in the stationary-distribution step, not in VFI itself.
- The error appears to come from a length mismatch inside a `sparse(...)` call in the semi-exogenous `nProbs` distribution code.

## Reproduction

In my application:

- finite horizon
- one endogenous asset state
- one decision variable `d`
- semi-exogenous states `semiz = (l, g)`
- no exogenous Markov state `z`
- permanent types handled through `*_PType` wrappers

The toolkit is called through:

- `ValueFnIter_Case1_FHorz_PType`
- `StationaryDist_Case1_FHorz_PType`

When I run with:

```matlab
vfoptions.gridinterplayer = 0;
```

the toolkit solve and stationary distribution complete successfully.

When I run with:

```matlab
vfoptions.gridinterplayer = 1;
vfoptions.ngridinterp = 15;
```

VFI completes, but the code fails during the stationary-distribution step.

## Error

The stack trace is:

```text
Error using sparse
Vectors must be the same length.

Error in StationaryDist_FHorz_Iteration_SemiExo_nProbs_noz_raw (line 47)
    Gammatranspose=sparse(Policy_aprimesemiz(:,:,jj),II2,PolicyProbs(:,:,jj),N_a*N_semiz,N_a*N_semiz);

Error in StationaryDist_FHorz_SemiExo (line 99)
    StationaryDist=StationaryDist_FHorz_Iteration_SemiExo_nProbs_noz_raw(...)

Error in StationaryDist_FHorz_Case1 (line 201)
    StationaryDist=StationaryDist_FHorz_SemiExo(...)

Error in StationaryDist_Case1_FHorz_PType (line 234)
    StationaryDist_ii=StationaryDist_FHorz_Case1(...)
```

Relevant toolkit functions:

- `StationaryDist_FHorz_SemiExo`
- `StationaryDist_FHorz_Iteration_SemiExo_nProbs_noz_raw`
- also likely the sibling `nProbs` semi-exogenous routines

## Why this seems to fail

Inside:

- `StationaryDist/FHorz/nProbs/StationaryDist_FHorz_Iteration_SemiExo_nProbs_noz_raw.m`

the code first compresses the semi-exogenous transition matrix to keep only the nonzero successors:

```matlab
N_semizshort = max(max(max(sum((pi_semiz_J>0),2))));
```

Then `Policy_aprimesemiz` and `PolicyProbs` are built using `N_probs * N_semizshort` entries per current state.

But `II2` is built using `N_probs * N_semiz`:

```matlab
II2 = repelem((1:1:N_a*N_semiz)',1,N_semiz*N_probs);
```

So if the semi-exogenous transition matrix is sparse and `N_semizshort < N_semiz`, then:

- `Policy_aprimesemiz(:,:,jj)` has width `N_probs * N_semizshort`
- `PolicyProbs(:,:,jj)` has width `N_probs * N_semizshort`
- `II2` has width `N_probs * N_semiz`

and the `sparse(...)` inputs no longer have matching lengths.

That seems consistent with the observed error.

## Proposed fix

In:

- `StationaryDist_FHorz_Iteration_SemiExo_nProbs_noz_raw.m`

replace:

```matlab
II2 = repelem((1:1:N_a*N_semiz)',1,N_semiz*N_probs);
```

with:

```matlab
II2 = repelem((1:1:N_a*N_semiz)',1,N_semizshort*N_probs);
```

## Likely related files

The same issue may also be present in the sibling semi-exogenous `nProbs` routines:

- `StationaryDist_FHorz_Iteration_SemiExo_nProbs_noz_e_raw.m`
- `StationaryDist_FHorz_Iteration_SemiExo_nProbs_e_raw.m`

Those files also appear to construct `II2` using `N_semiz * N_probs` even though the corresponding policy objects are expanded with `N_semizshort`.

For reference, the `z` version:

- `StationaryDist_FHorz_Iteration_SemiExo_nProbs_raw.m`

already uses `N_semizshort * N_probs` when building `II2`, which looks like the consistent implementation.

## Bottom line

My current understanding is:

- `gridinterplayer = 0` works
- `gridinterplayer = 1` fails only in the semi-exogenous stationary-distribution `nProbs` path
- the likely cause is an index-length mismatch caused by using `N_semiz` instead of `N_semizshort` when constructing `II2`

If useful, I can also provide a minimal reproducible example based on the life-cycle model that triggered this.
