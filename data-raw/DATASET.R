library(squat)
toto=spa:::learn_center_clust
# A=distDTW(toto, resample = FALSE, disable_normalization = TRUE)

qts11=toto[[1]]
qts11
qts12=toto[[1]]
qts2=toto[[2]]

bench::mark(
  a=derivative_qts(qts11)
)



bench::mark(
  A=distDTW(toto),
)

B=DTW(qts11, qts12)
B$normalizedDistance
C=DTW(qts11, qts2)
C$normalizedDistance
D=normalize_qts(qts11)
E=qts2angle(qts11)
G=reorient_qts(qts11)

x=rnorm(4)
x=rotations::as.Q4(rbind(x))
xx=as.numeric(x)
y=rnorm(4)
y=rotations::as.Q4(rbind(y))
yy=as.numeric(y)
xxx=xx[1:3] / sqrt(sum(xx[1:3]^2))
yyy=yy[1:3] / sqrt(sum(yy[1:3]^2))

bench::mark(
  new = inner_product_with_yinit(xx, yy),
  old=spa:::inner_product_with_yinit(xx, yy)
)

bench::mark(
  new = calibrate_xy(qts11, xx),
  old = spa:::calibrate_xy(qts11, xx)
)

bench::mark(
  new = rot_q(xxx, yyy),
  old = spa:::rot_q(xxx, yyy)
)

w=xx[1]
xx[1:3]=xx[2:4]
xx[4]=w
bench::mark(
  as.numeric(spa:::log.Q4(x)),
  log_quat(xx), check = FALSE
)


old = smooth_qts(qts11, alpha = 0.8)
new = smooth_qts(qts11)
spa:::plot_qts(list(qts11, old, new), clust = as.factor(1:3))
