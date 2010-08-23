#!/bin/sh


awk '(NR==1){print $1, -log($2), $3/$2}' avgent; var1=$(awk '($1==1){print -log($2)}' avgent) ; err1=$(awk '($1==1){print $3*$3/($2*$2)}' avgent) ;

awk '(NR==2 ){print $1, -log($2)+v1, sqrt($3*$3/($2*$2) + e1)}' v1=$var1 e1=$err1 01/avgent ; var1=$(awk '($1==2 ){print -log($2)+v1}' v1=$var1 01/avgent ) ; err1=$(awk '($1==2 ){print $3*$3 + e1}' e1=$err1 01/avgent ) ;
awk '(NR==3 ){print $1, -log($2)+v1, sqrt($3*$3/($2*$2) + e1)}' v1=$var1 e1=$err1 02/avgent ; var1=$(awk '($1==3 ){print -log($2)+v1}' v1=$var1 02/avgent ) ; err1=$(awk '($1==3 ){print $3*$3 + e1}' e1=$err1 02/avgent ) ;
awk '(NR==4 ){print $1, -log($2)+v1, sqrt($3*$3/($2*$2) + e1)}' v1=$var1 e1=$err1 03/avgent ; var1=$(awk '($1==4 ){print -log($2)+v1}' v1=$var1 03/avgent ) ; err1=$(awk '($1==4 ){print $3*$3 + e1}' e1=$err1 03/avgent ) ;
awk '(NR==5 ){print $1, -log($2)+v1, sqrt($3*$3/($2*$2) + e1)}' v1=$var1 e1=$err1 04/avgent ; var1=$(awk '($1==5 ){print -log($2)+v1}' v1=$var1 04/avgent ) ; err1=$(awk '($1==5 ){print $3*$3 + e1}' e1=$err1 04/avgent ) ;
awk '(NR==6 ){print $1, -log($2)+v1, sqrt($3*$3/($2*$2) + e1)}' v1=$var1 e1=$err1 05/avgent ; var1=$(awk '($1==6 ){print -log($2)+v1}' v1=$var1 05/avgent ) ; err1=$(awk '($1==6 ){print $3*$3 + e1}' e1=$err1 05/avgent ) ;
awk '(NR==7 ){print $1, -log($2)+v1, sqrt($3*$3/($2*$2) + e1)}' v1=$var1 e1=$err1 06/avgent ; var1=$(awk '($1==7 ){print -log($2)+v1}' v1=$var1 06/avgent ) ; err1=$(awk '($1==7 ){print $3*$3 + e1}' e1=$err1 06/avgent ) ;
awk '(NR==8 ){print $1, -log($2)+v1, sqrt($3*$3/($2*$2) + e1)}' v1=$var1 e1=$err1 07/avgent ; var1=$(awk '($1==8 ){print -log($2)+v1}' v1=$var1 07/avgent ) ; err1=$(awk '($1==8 ){print $3*$3 + e1}' e1=$err1 07/avgent ) ;
awk '(NR==9 ){print $1, -log($2)+v1, sqrt($3*$3/($2*$2) + e1)}' v1=$var1 e1=$err1 08/avgent ; var1=$(awk '($1==9 ){print -log($2)+v1}' v1=$var1 08/avgent ) ; err1=$(awk '($1==9 ){print $3*$3 + e1}' e1=$err1 08/avgent ) ;
awk '(NR==10){print $1, -log($2)+v1, sqrt($3*$3/($2*$2) + e1)}' v1=$var1 e1=$err1 09/avgent ; var1=$(awk '($1==10){print -log($2)+v1}' v1=$var1 09/avgent ) ; err1=$(awk '($1==10){print $3*$3 + e1}' e1=$err1 09/avgent ) ;
awk '(NR==11){print $1, -log($2)+v1, sqrt($3*$3/($2*$2) + e1)}' v1=$var1 e1=$err1 10/avgent ; var1=$(awk '($1==11){print -log($2)+v1}' v1=$var1 10/avgent ) ; err1=$(awk '($1==11){print $3*$3 + e1}' e1=$err1 10/avgent ) ;
awk '(NR==12){print $1, -log($2)+v1, sqrt($3*$3/($2*$2) + e1)}' v1=$var1 e1=$err1 11/avgent ; var1=$(awk '($1==12){print -log($2)+v1}' v1=$var1 11/avgent ) ; err1=$(awk '($1==12){print $3*$3 + e1}' e1=$err1 11/avgent ) ;
awk '(NR==13){print $1, -log($2)+v1, sqrt($3*$3/($2*$2) + e1)}' v1=$var1 e1=$err1 12/avgent ; var1=$(awk '($1==13){print -log($2)+v1}' v1=$var1 12/avgent ) ; err1=$(awk '($1==13){print $3*$3 + e1}' e1=$err1 12/avgent ) ;
awk '(NR==14){print $1, -log($2)+v1, sqrt($3*$3/($2*$2) + e1)}' v1=$var1 e1=$err1 13/avgent ; var1=$(awk '($1==14){print -log($2)+v1}' v1=$var1 13/avgent ) ; err1=$(awk '($1==14){print $3*$3 + e1}' e1=$err1 13/avgent ) ;
awk '(NR==15){print $1, -log($2)+v1, sqrt($3*$3/($2*$2) + e1)}' v1=$var1 e1=$err1 14/avgent ; var1=$(awk '($1==15){print -log($2)+v1}' v1=$var1 14/avgent ) ; err1=$(awk '($1==15){print $3*$3 + e1}' e1=$err1 14/avgent ) ;
awk '(NR==16){print $1, -log($2)+v1, sqrt($3*$3/($2*$2) + e1)}' v1=$var1 e1=$err1 15/avgent ; var1=$(awk '($1==16){print -log($2)+v1}' v1=$var1 15/avgent ) ; err1=$(awk '($1==16){print $3*$3 + e1}' e1=$err1 15/avgent ) ;
awk '(NR==17){print $1, -log($2)+v1, sqrt($3*$3/($2*$2) + e1)}' v1=$var1 e1=$err1 16/avgent ; var1=$(awk '($1==17){print -log($2)+v1}' v1=$var1 16/avgent ) ; err1=$(awk '($1==17){print $3*$3 + e1}' e1=$err1 16/avgent ) ;
awk '(NR==18){print $1, -log($2)+v1, sqrt($3*$3/($2*$2) + e1)}' v1=$var1 e1=$err1 17/avgent ; var1=$(awk '($1==18){print -log($2)+v1}' v1=$var1 17/avgent ) ; err1=$(awk '($1==18){print $3*$3 + e1}' e1=$err1 17/avgent ) ;
awk '(NR==19){print $1, -log($2)+v1, sqrt($3*$3/($2*$2) + e1)}' v1=$var1 e1=$err1 18/avgent ; var1=$(awk '($1==19){print -log($2)+v1}' v1=$var1 18/avgent ) ; err1=$(awk '($1==19){print $3*$3 + e1}' e1=$err1 18/avgent ) ;
awk '(NR>19 ){print $1, -log($2)+v1, sqrt($3*$3/($2*$2) + e1)}' v1=$var1 e1=$err1 19/avgent ;

