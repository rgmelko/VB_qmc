#!/bin/sh


awk '(NR==1){print $1, -log($2), $3/$2}' avgent; var1=$(awk '($1==1){print -log($2)}' avgent) ; err1=$(awk '(NR==1){print $3*$3/($2*$2)}' avgent) ;

awk '(NR==2 ){print $1, -log($2)+v1, sqrt($3*$3/($2*$2) + e1)}' v1=$var1 e1=$err1 01/avgent ; var1=$(awk '(NR==2 ){print -log($2)+v1}' v1=$var1 01/avgent ) ; err1=$(awk '($1==2 ){print $3*$3/($2*$2) + e1}' e1=$err1 01/avgent ) ;
awk '(NR==3 ){print $1, -log($2)+v1, sqrt($3*$3/($2*$2) + e1)}' v1=$var1 e1=$err1 02/avgent ; var1=$(awk '(NR==3 ){print -log($2)+v1}' v1=$var1 02/avgent ) ; err1=$(awk '($1==3 ){print $3*$3/($2*$2) + e1}' e1=$err1 02/avgent ) ;
awk '(NR==4 ){print $1, -log($2)+v1, sqrt($3*$3/($2*$2) + e1)}' v1=$var1 e1=$err1 03/avgent ; var1=$(awk '(NR==4 ){print -log($2)+v1}' v1=$var1 03/avgent ) ; err1=$(awk '($1==4 ){print $3*$3/($2*$2) + e1}' e1=$err1 03/avgent ) ;
awk '(NR==5 ){print $1, -log($2)+v1, sqrt($3*$3/($2*$2) + e1)}' v1=$var1 e1=$err1 04/avgent ; var1=$(awk '(NR==5 ){print -log($2)+v1}' v1=$var1 04/avgent ) ; err1=$(awk '($1==5 ){print $3*$3/($2*$2) + e1}' e1=$err1 04/avgent ) ;
awk '(NR>5  ){print $1, -log($2)+v1, sqrt($3*$3/($2*$2) + e1)}' v1=$var1 e1=$err1 05/avgent ;

