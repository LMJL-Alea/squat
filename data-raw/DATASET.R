toto=spa:::learn_center_clust
qts11=toto[[1]]
qts12=toto[[1]]
qts2=toto[[2]]
A=distDTW(toto)
B=DTW(qts11, qts12)
B$normalizedDistance
C=DTW(qts11, qts2)
C$normalizedDistance
D=normalize_qts(qts11)
E=qts2angle(qts11)
G=reorient_qts(qts11)
