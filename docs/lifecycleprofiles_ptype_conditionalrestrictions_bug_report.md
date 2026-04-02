# Possible bug in `LifeCycleProfiles_FHorz_Case1_PType`

I think there is a bug in

`EvaluateFnOnAgentDist/FHorz/PType/LifeCycleProfiles_FHorz_Case1_PType.m`

when all of these are used together:

- `PType`
- `simoptions.conditionalrestrictions`
- `simoptions.groupptypesforstats = 1`

I get:

```text
Error using accumarray
Second input must be a vector with one element for each row in the first input, or a scalar.

Error in LifeCycleProfiles_FHorz_Case1_PType (line 1117)
AllRestrictedWeights_rrffjj=accumarray(sortindex,AllRestrictedWeights_rrffjj,[],@sum);
```

The issue seems to be this:

Earlier, in the grouped unconditional-stats block, the code does

```matlab
[AllValues.(FnsToEvalNames{ff}).(jgroupstr{jj}),~,sortindex]=unique(AllValues.(FnsToEvalNames{ff}).(jgroupstr{jj}));
AllWeights.(FnsToEvalNames{ff}).(jgroupstr{jj})=accumarray(sortindex,AllWeights.(FnsToEvalNames{ff}).(jgroupstr{jj}),[],@sum);
```

so `AllValues.(...)` gets overwritten by its unique support points.

Later, in the grouped conditional-restriction block, the code does

```matlab
AllValues_rrffjj=gather(AllValues.(FnsToEvalNames{ff}).(jgroupstr{jj})(:));
AllRestrictedWeights_rrffjj=gather(AllRestrictedWeights.(CondlRestnFnNames{rr}).(FnsToEvalNames{ff}).(jgroupstr{jj})(:));
[AllValues_rrffjj,~,sortindex]=unique(AllValues_rrffjj);
AllRestrictedWeights_rrffjj=accumarray(sortindex,AllRestrictedWeights_rrffjj,[],@sum);
```

But at that point:

- `AllValues.(...)` has already been collapsed by `unique(...)`
- `AllRestrictedWeights.(...)` is still at the old pre-collapsed length

So `sortindex` and `AllRestrictedWeights_rrffjj` no longer have matching lengths, which seems to explain the `accumarray` error.

My guess is that the grouped unconditional-statistics block should use temporary variables rather than overwrite `AllValues.(...)`, or else the grouped conditional-restriction block needs to rebuild matching values/weights from the same pre-collapsed objects.

This looks like a generic bookkeeping bug in the grouped-`PType` conditional-restriction path, not something model-specific.
