#########library "plsRbeta"
library(plsRglm)
library(plsRbeta)
library(tidyverse)

##  
## plsRbetaの説明が不親切なので勉強する
## ベータ回帰の場合、予測や誤差の評価が解析的計算できない？らしく、基本的に
## ブーツストラップを使う。ここではパッケージ公式のガソリンのチュートリアルから
## 勉強してみる
##

data("GasolineYield",package="betareg")
yGasolineYield <- GasolineYield$yield
XGasolineYield <- GasolineYield[,2:5]
hist(yGasolineYield)
plot(GasolineYield)

# running plsRbeta, link function is set at logit as default
# X is standardised as default
# nt specifies number of component to extract, which can be determined using AIC/BIC

# first, let's check full model with Cross validation
modpls_cv <- PLS_beta_kfoldcv_formula(yield~.,data=GasolineYield,nt=5,modele="pls-beta",verbose=FALSE)
full_cv <- kfolds2CVinfos_beta(modpls_cv) # extract output

#元論文によると、どの指標が適切かどうか定説はないらしい。
#AICとBICは大目に主成分を選ぶ傾向があるらしい。というわけで、表から組み合わせてこんなもんだろう、と選ぶのがいいらしい。
full_cv

##ここでkfolds2CVinfos_beta()で帰ってくる指標が何なのかわからなかった。
##githubでコードも見てみたけど、解読できず。
cv_results <- modpls_cv$results_kfolds # list of predicted value (pred_Y) using different numbers of comp (column) 
cv_data <- modpls_cv$dataY_kfolds # obs_Y. Response variable. Exactly same as 'yGasolineYield'

df_results <- cv_results %>% 
  unlist() %>% 
  matrix(nrow=5,ncol=32) %>% 
  t() %>% 
  as.data.frame()

sum((yGasolineYield-df_results$V2)^2/var(yGasolineYield))

##　基本的なモデルは普通のformulaで書ける。モデルをpls-betaと指定するとデフォルトでlogitがリンク関数になる
## Google Earth Engineで予測したいときはモデル$Coeffsで標準化する前のオリジナルのスケールでの係数が得られる。
#running pls beta
modpls <- plsRbeta(yield~., data=GasolineYield, nt=3, modele="pls-beta", verbose=FALSE)
modpls$Coeffs # likely coefficients used with the original data scale
#modpls$PredictY[1,] # standardized X by subtracting scaled:center then divided by scaled:scale

## ブーツストラップはbootパッケージを使用しているらしい。typeboot="fmodel_np"は高速かつ安定した出力が得られるとのこと。
GazYield.boot <- bootplsbeta(modpls, typeboot="fmodel_np", sim="ordinary", stype="i", R=250)
GazYield.boot # originalが平均値をbias補正したもの

# boxplotでブーツストラップの分布を図示できる
plsRglm::boxplots.bootpls(GazYield.boot)
plsRglm::plots.confints.bootpls(plsRglm::confints.bootpls(GazYield.boot))

# 信頼区間の数値を得る（=パラメトリックに分布を近似またはノンパラな分位点）にはいくつかの方法があるが
# BCaというのが無難そう。これはconfints.bootpls関数の７、８列目で得られる。
temp.ci <- confints.bootpls(GazYield.boot,indice=1:4)
plots.confints.bootpls(temp.ci,typeIC="BCa", legendpos ="topright") # BCa (Bias Corrected and acceralated) should be fine in most cases
ind.BCa.GazYT1 <- (temp.ci[,7]<0&temp.ci[,8]<0)|(temp.ci[,7]>0&temp.ci[,8]>0) # Normal Lower then Upper Limit, Basic Lower then Upper Limit, Percentile Lower then Upper Limit, BCa Lower then Upper Limit.
ind.BCa.GazYT1

## ループで異なる主成分数のbootstrapを回せばどの変数が選択されるかわかる。
## Then check with other numbers of COMPs
output_ind_BCa <- data.frame(matrix(rep(NA, 4), nrow=1))[numeric(0), ]
colnames(output_ind_BCa) <- names(ind.BCa.GazYT1[1:4])

for(i in 1:5){
  temp_plsbeta <- plsRbeta(yield~., data=GasolineYield, nt=i, modele="pls-beta", verbose=FALSE)
  temp_boot <- bootplsbeta(temp_plsbeta,typeboot="fmodel_np",R=1000)
  temp_ci <- confints.bootpls(temp_boot,indices=1:4)
  temp_ind <- (temp.ci[,7]<0&temp.ci[,8]<0)|(temp.ci[,7]>0&temp.ci[,8]>0)
  output_ind_BCa <- rbind(output_ind_BCa, temp_ind)
}

output_ind_BCa
