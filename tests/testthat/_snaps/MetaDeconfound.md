# logistic regression

    Code
      MetaDeconfound(featureMat = feature[, c("MS0035"), drop = F], metaMat = metaMat,
      logLevel = "WARN", returnLong = T, logistic = T)
    Output
      [] 
      	Unallowed characters detected in rownames and/or colnames of featureMat and/or metaMat!
      	metadeconfoundR will try to remove these characters using the make.names() function.
      [] Separation for: MS0035 and Status
      [] Separation for: MS0035 and Dataset
      [] Separation for: MS0035 and Metformin
        feature     metaVariable        Ps        Qs          Ds status
      1  MS0035           Status        NA        NA  0.03130230     NS
      2  MS0035          Dataset        NA        NA         Inf     NS
      3  MS0035        Metformin        NA        NA  0.01971322     NS
      4  MS0035 continuous_dummy 0.4091398 0.4091398 -0.34573003     NS
      5  MS0035    altered_dummy 0.6115777 0.6115777 -0.24104683     NS

