# Possible bug in `EvalFnOnAgentDist_ValuesOnGrid_FHorz_Case1_PType` for finite-horizon `PType` models with `semiz` and no `z`

I think there is a bug in

`EvaluateFnOnAgentDist/FHorz/PType/EvalFnOnAgentDist_ValuesOnGrid_FHorz_Case1_PType.m`

when it is used in a finite-horizon `PType` model with:

- no ordinary exogenous state (`n_z = 0`)
- a semi-exogenous state `semiz`

## Symptom

Calling

```matlab
ValuesOnGrid = EvalFnOnAgentDist_ValuesOnGrid_FHorz_Case1_PType( ...
    Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,N_i,d_grid,a_grid,z_grid,simoptions);
```

fails with

```text
Error using gpuArray/reshape
Number of elements must not change.

Error in EvalFnOnAgentDist_ValuesOnGrid_FHorz_Case1_PType (line 253)
ValuesOnGrid.(FnNames{kk}).(Names_i{ii})=gather(reshape(ValuesOnGrid_ii.(FnNames{kk}),[n_a_temp,N_j_temp]));
```

## Why this looks like a bug

In this setup, `n_z = 0`, but `simoptions.n_semiz` is nonzero. So the evaluated objects should still live on the grid

```text
(a, semiz, j)
```

not just

```text
(a, j)
```

The non-`PType` helper

`EvaluateFnOnAgentDist/FHorz/EvalFnOnAgentDist_ValuesOnGrid_FHorz_Case1.m`

handles this by folding `semiz` into the effective exogenous grid before reshaping. So with `semiz`, the output should have size like

```text
[n_a_temp, n_semiz_temp, N_j_temp]
```

or, in the more general notation used in the wrapper,

```text
[n_a_temp, n_ze_temp, N_j_temp]
```

However, in the `PType` wrapper, the code falls into the branch

```matlab
if prod(n_ze_temp)==0
```

and then reshapes to

```matlab
[n_a_temp, N_j_temp]
```

which is only valid when there are truly no exogenous or semi-exogenous states.

But with `semiz`, `ValuesOnGrid_ii.(FnNames{kk})` contains values over `(a, semiz, j)`, so it has more than `n_a_temp * N_j_temp` elements. That is why the reshape fails.

## Suspected issue

It looks like the `PType` wrapper is incorrectly treating the case

- `n_z = 0`
- `n_semiz > 0`

as if there were no exogenous-state dimensions at all.

## Likely fix

The wrapper should preserve the `semiz` dimension in the same way as the non-`PType` helper does, so that this case reshapes to something like

```matlab
[n_a_temp, n_ze_temp, N_j_temp]
```

rather than

```matlab
[n_a_temp, N_j_temp]
```

when `semiz` is present.

