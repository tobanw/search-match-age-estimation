# Estimation notes

* TODO: if pooling the lower boundary population, need to adjust aging rate to reflect that only 1/k of the initial age group will age out of it
	* Currently, not pooling at all; could just acknowledge the boundary bias problem but note that it doesn't affect the action in 30s and onward
	* for edu, just start model at 24 with inflows from 23, i.e., `min.age = 23` in SQL queries
* GMM estimation drops terminal age moments, so meeting rates can't be biased by the edge cases
* moment fitting insights: MF (DF) strictly monotonic in lambda (delta)
	* alpha decreases in lambda, but MF still increases in lambda
	* optimal delta increases in alpha, so optimal xi and delta are inversely related
	* MF increases a lot in delta, hardly in lambda
* model fit: ***can only fit MF/DF values that are consistent with steady states for u,m!***
	* can't fit huge MF for young: increasing lambda is largely offset by decreased alpha as stock of marriages too small
	* probably a result of migration outflows for couples: the observed MF can only be accounted for by excessively large DF, unless there are migration outflows of couples
* dropped edu as it isn't stationary over cohorts; but then loss shot up 30x!
	* likely the squared difference terms: fewer but larger differences squared
	* also worse fit, can't capture the differing marriage patterns by type
* type specification
	* male marriageability: binary type based on an index of income, employment, college

## Non-parametric estimates

* age-only model: why is xi so large? makes alpha tiny and hence s is a large negative and dominates prod
* prod wackiness: from alpha bottoming out due to truncation, try adjusting truncation threshold
	* steep cliffs in alpha come purely from differencing out the marriage inflows (due to marriage density spreading out in age-gap)
	* these cliffs or sidewalls show up when mar.bw is around 3
* prod choppy: first diffs? f = (s - d1 s1) + (v - dd1 v1) - delta d c1 mu(alpha)
	* smoking gun: (s - d1 s1) has same pattern as prod! As s is just inverse of alpha, this implicates double diff of marriage stocks
	* this term dominates when alpha is tiny (s is a large negative), which occurs with large xi
	* also explains 20s prod pit: alpha/surplus is rapidly rising and peaks in 30s, so s diff is negative
* value functions rise until mid-20s, then decline uniformly
	* expect them to monotonically decrease; this is likely a boundary bias issue because not including the mass of younger people available in the marriage market

## Smoothing: non-parametric local linear regression

### Custom uncontrained bandwidth kernel

* local-poly takes 5 min for linear, 7 min for quadratic, 12 min for cubic
* (5, 4) for all three (pop, flow, mig): needs more smoothing
	* DF is jumpy and pretty off (peaks too late)
	* prod boiling, alpha could be a bit smoother
	* mar stock and flow good, mig pretty rough
* (6, 5.5) for all three (pop, flow, mig): needs wider bw span, e.g. 6-5
	* DF is jumpy and pretty off
	* prod boiling, alpha worse
	* mar stock and flow good, mig is visibly picking up diagonal striations
* (10, 8) for all three (pop, flow, mig): better ratio?
	* DF miss doesn't seem as bad
	* prod/alpha not much better than 8.5
* (10, 8.5) for all three (pop, flow, mig): better ratio?
	* prod smoother, alpha good
* (10, 9) for all three (pop, flow, mig): better ratio
	* DF smoother but still pretty off
	* prod choppy but can see overall pattern, alpha better
	* mar stock and flow good, mig is better
* (10, 9.5) for all three (pop, flow, mig)
	* DF smoother but still pretty off
	* prod choppy but can see overall pattern, alpha not great
	* mar stock and flow good, mig is visibly picking up diagonal striations
* (16, 13) for all three (pop, flow, mig)
	* prod much better, alpha good
	* mig looks good

### Age-only model (18-65)

* individual objects are nice and smooth with CV bandwidth
* CV bandwidth for marriage masses is 0.75-0.8 across MSAs
	* too rough, need heavier smoothing for marriage objects
* xi=0.82; this makes alpha very choppy, resulting in garbled production and crazy spiky DF
* bw=1.5: xi jumps up to 2.17, DF settles down (pretty nice fit), but production is still garbled
* bw=2.5: xi jumps up to 26.5, alpha looks pretty nice (large xi makes alpha very small), but production still sucks (get an archipelago of islands dotting the a=b line); not sure why, as everything else looks nice and smooth (pops, mig, val, alpha)
* testing different bw on just the top 5 msa's: xi is much higher (suggesting IRS tech?)
* bw=1.6: alpha has the diagonal islands, prod is roiling bubbles, pops are smooth except for marriage outflows which have islands
* bw=2: prod is smoother but not enough, alpha and marriage outflows still have islands
* bw=3: prod is now smooth with lumps, everything is smooth and outflows are only mildly lumpy
* try different smoothing for marriage mass/flows/outflows (mass and outflows might need most because of double diff)
	* Notes: both arrival rates are increasing in bw, NYC tends to overpredict MF
	* higher bw makes prod smoother/wider but doesn't fundamentally change existence of 20s pit, brings sidewall artifacts
	* MF fit: stock/flow bws should be calibrated so that MF spread matches up with model (lambda u alpha); e.g., if bws too different, then moment won't match because of spread
	* bw=2 for each: model MF too sharp, overpredicts along a=b diagonal; need higher mar.bw so that model MF is more spread out
	* 3, 1.6, 2.2; alpha pretty smooth but prod is lumpy, mig is okay
	* 3, 1.6, 3; MF fit bad - 3 vs 1.6 too large of a gap, alpha looks smoother but prod not much better, mig is better
	* ageonly4: 4, 1.6, 4; model MF is too spread (mar.bw too large relative to flow.bw), the sidewall artifacts are appearing, but prod is a lot smoother
	* ageonly5: 5, 1.6, 5; low prod for young shows up
	* 2.5, 2, 3; MF fit better but still overestimating
	* 2.5, 1.5, 2.5: MF fit underestimates a sliver along a=b, but overestimates more generally around it
	* 2.5, 1.7, 2.5: MF fit similar though not as bad in the sliver
	* 3, 2.5, 3; MF fit still overestimating, but might be best so far (150 max colorbar, no sliver)
	* 2.2, 1.8, 2.2; MF fit overestimating but no sliver
	* ageonly2: 2, 1.6, 2; MF fit overestimates even more
* ageonly25, on top 20: 2.5, 2, 2.5; prod roils too much, need heavier smoothing
* ageonly35: 3.5, 3, 3.5; prod much nicer! but MF fit no better
* ageonly4: 4, 3.5, 4; prod even more sensible! MF fit a bit better
* ageonly5: 5, 4.5, 5; ?
* ageonly3, on top 3: 3, 2.5, 3; doesn't help much with MF fit, so it's not the NYC vs Dallas tension
* ageonly, on top 20: 3, 2.5, 3; top 20 vs top 3 doesn't affect rates much (75 vs 81) nor other objects

### With edu/race types (25-65)

* bw=2 works well for thick (x,y) pairs, but leaves the rest lumpy
* bw=3 works well for thin (x,y) pairs too
* even bw=3 gives garbled f... try heavier smoothing for singles too
* trying bw=6...

## Assumptions and their implications

### Steady state

* cities are pretty dynamic and change over decades, so age structure might be way off steady state
* also general social changes over time in marriage/family
* college attendance has risen a lot, especially for women; this shows up very clearly in the total populations of men/women by type

### Continuous time model but annual data

* assuming infinitesimal time interval (dt) and hence no simultaneous events
	* but clearly aging and other events must occur together in data interval
* arrival rates are defined per unit of time (annual), but multiplied by dt -> 0 to use Ito and drop simultaneous events
* mainly a problem for flow balance equation (alpha), but maybe not as bad as it seems:
	* inflow is entire stock of younger couples (no subtraction of death/divorce)
	* outflow is entire stock of couples plus death/divorce: impossible
	* net aging flow is just difference of stocks, so the extra inflow is offset by the extra outflow of death/divorce; so net age flow is slighly biased and should be multiplied by (rho - death - divorce) and not just rho
	* could try subtracting death/divorce from aging rate to get true age flow, but likely a very small effect since death/divorce rates are small and move slowly

### Dating is instantaneous

* takes time to date before getting married, so sex ratio effects will show up in MF with a lag


## Optimization

* Monotonicity of moments (mostly) means no local minima, just use gradient-based optimization

### Flexible specification of arrival rates

* Provide gradient function? Compute pt-wise and then `vec()` it
* different xi's only interact with one another through delta: breaks curse of dim in optimization!
	* Alternate until convergence: 1) delta update step; 2) pt-wise update each xi
* however, xi and delta both affect alpha, so this may converge to alternating min (analogous to best responses sans randomization), but averaging them would give a good initial guess
	* Given delta, supply gradient for xi -- pointwise this will be the sum of the MF and DF parts
	* Or just solve grad = 0 analytically and provide current minimum in one shot

## Arrival rates

* identification of xi? objective function is extremely flat in xi, because model MF is proportional to `lambda alpha`, but lambda effectively cancels out when `delta m` is tiny relative to `lambda u`
	* lambda only has leverage when `delta m` is substantial: `lambda alpha = top / (u + delta/lambda m)`
	* `lambda = xi/sqrt(U)` is meeting rate: fraction of potential matches realized per unit of time; close to zero because U very large
	* because of scale of u (number of potential meetings), `lambda u` term can really blow up with lambda
	* `lambda u = xi sqrt(U) u/U `: meeting flow is fraction of type-pair potential meetings out of all potential meetings (`u/U`) times matching tech (`xi sqrt(U)`)
	* TESTING overweighting MF vs DF (equal weighting gives infinite xi as solution, because marginal gains to improving DF loss dominates)
		* tradeoff: overestimate MF or predict rise of divorce too late
		* I.e., lambda too high (overestimate MF) or too low (underestimate DF)
		* Implies model is not consistent with data: marriage and divorce are NOT governed by a common alpha!
		* Ageonly: equal = 7, 2MF = 1.37, 4MF = 0.6
* comparisons:
	* GJR: delta = 0.038, xi=0.1515
	* Gousse: delta = (0.025, 0.058), lambda (quadratic) = (0.0029, 0.0078)
