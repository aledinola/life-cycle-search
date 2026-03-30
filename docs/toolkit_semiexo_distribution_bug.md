# Possible Bug in `StationaryDist/FHorz/SemiExo/StationaryDist_FHorz_Iteration_SemiExo_noz_raw.m`

I think there is a bug in the way this function indexes the finite-horizon semi-exogenous transition matrix.

## Relevant object

The input `pi_semiz_J` is a 4-dimensional array with size

```matlab
[N_semiz, N_semiz, N_dsemiz, N_j]
```

where the dimensions are:

1. current semi-exogenous state
2. next semi-exogenous state
3. decision affecting the semi-exogenous transition
4. age

So in a finite-horizon model, the transition matrix is allowed to vary with age.

## Relevant code

In `StationaryDist_FHorz_Iteration_SemiExo_noz_raw.m`:

```matlab
Policy_dsemiexo=reshape(Policy_dsemiexo,[N_a*N_semiz,1,N_j]);

semizindex=repelem(gpuArray(1:1:N_semiz)',N_a,1) ...
    + N_semiz*gpuArray(0:1:N_semiz-1) ...
    + (N_semiz*N_semiz)*(Policy_dsemiexo-1);
```

and later:

```matlab
semiztransitions=gather(pi_semiz_J(semizindex(:,:,jj)));
```

## Why this looks wrong

The array `semizindex` only encodes linear indexing for the first 3 dimensions of `pi_semiz_J`:

- current semiz
- next semiz
- decision

It does **not** include any offset for the 4th dimension, age.

In particular, there is no term of the form

```matlab
(N_semiz*N_semiz*N_dsemiz)*(jj-1)
```

in the construction of `semizindex`.

Therefore, `semizindex(:,:,jj)` is an index for a 3-dimensional array of size

```matlab
[N_semiz, N_semiz, N_dsemiz]
```

but it is then being used to index the 4-dimensional array `pi_semiz_J`.

## Why this implies a bug

If `pi_semiz_J` varies with age, then the code should select the age-`jj` slice of the transition matrix before applying the linear index over `(semiz, semiz', d)`.

But the current code does:

```matlab
pi_semiz_J(semizindex(:,:,jj))
```

instead of first restricting to age `jj`.

So the age dimension is not being indexed correctly.

## What I think the code should do instead

It seems the safe way is:

```matlab
pi_semiz_jj = pi_semiz_J(:,:,:,jj);
semiztransitions = gather(pi_semiz_jj(semizindex(:,:,jj)));
```

because `semizindex(:,:,jj)` is exactly an index for the 3-dimensional object

```matlab
pi_semiz_jj(semiz, semizprime, dsemiz).
```

## Bottom line

So the issue is:

- `pi_semiz_J` is 4-dimensional
- `semizindex` only indexes 3 dimensions
- the current code applies the 3-dimensional linear index directly to the 4-dimensional array

This appears to ignore the finite-horizon age dimension of the semi-exogenous transition matrix.
