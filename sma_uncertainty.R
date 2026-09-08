## AUTO-GENERATED uncertainty inputs (DSA ranges + PSA distributions)
## Source: SMA Excel Model.xlsm (Utility, Costs, TP sheets). Regenerate; do not hand-edit.

## Utilities ~ Beta(alpha,beta)
util_beta <- list(
  NS=c(a=93.00843749999997, b=227.7103125),
  S=c(a=25.311108033241002, b=39.589168975069256),
  W=c(a=32.19428571428571, b=12.520000000000003),
  PV=c(a=81.03500000000001, b=345.46500000000003)
)

## State health-state cost coefficient of variation (SE/mean); PSA ~ Gamma
cost_cv <- c(NS=0.2433, S=0.2998, W=0.1167, PV=0.2433)
## Drug-cost CV per arm (PSA ~ Gamma; DSA +/-20%)
drug_cv <- c(Nusinersen=0.2, OA=0.05, Risdiplam=0.20)  ## risdiplam CV: placeholder, set from source

## Fixed transition PSA/DSA: NS-row probs ~ Beta(method-of-moments from base,SE);
##   S_W (OA only) ~ Beta(a,b). SE derived from 95% CI (Lo/Up) or SE column.
tp_unc <- list(
  BSC = list(NS_S=list(kind="prob", mean=0.0, se=0.00263146), NS_PV=list(kind="prob", mean=0.16308168, se=0.16419083), NS_Death=list(kind="prob", mean=0.12944895, se=0.15618259)),
  Nusinersen = list(NS_S=list(kind="prob", mean=0.02849779, se=0.02156061), NS_PV=list(kind="prob", mean=0.0913675, se=0.08992538), NS_Death=list(kind="prob", mean=0.04920689, se=0.09071578)),
  OA = list(NS_S=list(kind="prob", mean=0.12998521, se=0.08281789), NS_PV=list(kind="prob", mean=0.00971198, se=0.0004856), NS_Death=list(kind="prob", mean=0.00971198, se=0.0004856), S_W=list(kind="beta", a=389.1755, b=14028.1006)),
  Risdiplam = list(NS_S=list(kind="prob", mean=0.10371905, se=0.0302213), NS_PV=list(kind="prob", mean=0.01552434, se=0.00077622), NS_Death=list(kind="prob", mean=0.0130834, se=0.00065417))
)
