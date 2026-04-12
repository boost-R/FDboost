# \`FDboost_fac\` S3 class for factorized FDboost model components

Model factorization with \`factorize()\` decomposes an \`FDboost\` model
into two objects of class \`FDboost_fac\` - one for the response and one
for the covariate predictor. The first is essentially an \`FDboost\`
object and the second an \`mboost\` object, however, in a 'read-only'
mode and slightly adjusted methods (method defaults).

## See also

\[factorize(), factorize.FDboost()\]
