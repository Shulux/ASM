###### Practical 2
###### Mushrooms

# We create a matrix with the data in.
mushroom_data <- matrix( c(101, 399, 57, 487,
                           12, 389, 150, 428 ), byrow = TRUE, ncol = 4 )

# Add dimension names as follows.
dimnames( mushroom_data ) <- list( Edibility = c("Edible", "Poisonous"),
                                   Cap_Shape = c("bell", "flat", "knobbed",
                                                 "convex/conical") )
# Have a look.
mushroom_data


# Let's look at the proportion of each shape of mushroom that are
# edible or poisonous.
mushroom_table <- addmargins( mushroom_data, 2, FUN = mean )
mushroom_table <- prop.table( mushroom_table, 2 )
mushroom_table <- addmargins(mushroom_table, 1 )
mushroom_table

# Let's perform a chi-square test.
chisq.test( mushroom_data )

# Barplot
barplot( prop.table( mushroom_data, margin = 2 ), density = 50,
         main = "Comparison of Edibility by Cap Shape", cex.main = 0.8,
         xlab = "treatment outcome", ylab = "proportions",
         legend.text = T, args.legend = list( x = 5.1, y = 1.25 ) )


# Sieve diagram
sieve( mushroom_data )


# Much can be drawn form this diagram. For example,
# we can see that there are far more edible bell mushrooms in our
# sample than would be expected under an assumption of independence.


chisq_Mush <- chisq.test( mushroom_data )
chisq_Mush

chisq_Mush$residuals
chisq_Mush$stdres

G2 <- function( data ){
  # computes the G2 test of independence
  # for a two-way contingency table of
  # data: IxJ matrix
  X2 <- chisq.test( data )
  Ehat <- X2$expected
  df <- X2$parameter
  term.G2 <- data * log( data / Ehat )
  term.G2[data==0] <- 0 # Because if data == 0, we get NaN
  Gij <- 2 * term.G2 # Individual cell contributions to G2 statistic.
  dev_res <- sign( data - Ehat ) * sqrt( abs( Gij ) )
  G2 <- sum( Gij ) # G2 statistic
  p <- 1 - pchisq( G2, df )
  return( list( G2 = G2, df = df, p.value = p,
                Gij = Gij, dev_res = dev_res ) )
}


G2mush <- G2( mushroom_data )
G2mush$p.value
G2mush$dev_res
G2mush$Gij

OR <- function( M ){ ( M[1,1] * M[2,2] ) / ( M[1,2] * M[2,1] ) }


nominal_OR <- function( M, ref_x = nrow( M ), ref_y = ncol( M ) ){
  # I and J
  I <- nrow(M)
  J <- ncol(M)
  # Odds ratio matrix.
  OR_reference_IJ <- matrix( NA, nrow = I, ncol = J )
  for( i in 1:I ){
    for( j in 1:J ){
      OR_reference_IJ[i,j] <- OR( M = M[c(i,ref_x), c(j,ref_y)] )
    }
  }
  OR_reference <- OR_reference_IJ[-ref_x, -ref_y, drop = FALSE]
  return(OR_reference)
}

nominal_OR(mushroom_data)

mosaic( mushroom_data,
        gp = shading_hcl,
        residuals_type = "Pearson" )


linear.trend <- function( table, x, y ){
  # linear trend test for a 2-way table
  # PARAMETERS:
  # freq: vector of the frequencies, given by rows
  # NI: number of rows
  # NJ: number of columns
  # x: vector of row scores
  # y: vector of column scores
  # RETURNS:
  # r: Pearson’s sample correlation
  # M2: test statistic
  # p.value: two-sided p-value of the asymptotic M2-test
  NI <- nrow( table )
  NJ <- ncol( table )
  rowmarg <- addmargins( table )[,NJ+1][1:NI]
  colmarg <- addmargins( table )[NI+1,][1:NJ]
  n <- addmargins( table )[NI+1,NJ+1]
  xmean <- sum( rowmarg * x ) / n
  ymean <- sum( colmarg * y ) / n
  xsq <- sqrt( sum( rowmarg * ( x - xmean )^2 ) )
  ysq <- sqrt( sum( colmarg * ( y - ymean )^2 ) )
  r <- sum( ( x - xmean ) %*% table %*% ( y - ymean ) ) / ( xsq * ysq )
  M2 = (n-1)*r^2
  p.value <- 1 - pchisq( M2, 1 )
  return( list( r = r, M2 = M2, p.value = p.value ) )
}


DoseResult <- matrix( c(47, 25, 12,
                        36, 22, 18,
                        41, 60, 55), byrow = TRUE, ncol = 3 )
dimnames( DoseResult ) <- list( Dose=c("High",
                                       "Medium",
                                       "Low"),
                                Result=c("Success",
                                         "Partial",
                                         "Failure") )
DoseResultTable <- addmargins(DoseResult)


local_OR <- function( M ){
  # I and J
  I <- nrow(M)
  J <- ncol(M)
  # Odds ratio matrix.
  OR_local <- matrix( NA, nrow = I-1, ncol = J-1 )
  for( i in 1:(I-1) ){
    for( j in 1:(J-1) ){
      OR_local[i,j] <- OR( M = M[c(i,i+1), c(j,j+1)] )
    }
  }
  return(OR_local)
}

global_OR <- function( M ){
  # I and J
  I <- nrow(M)
  J <- ncol(M)
  # Odds ratio matrix.
  OR_global <- matrix( NA, nrow = I-1, ncol = J-1 )
  for( i in 1:(I-1) ){
    for( j in 1:(J-1) ){
      OR_global[i,j] <- OR( M = matrix( c( sum( M[1:i,1:j] ),
                                           sum( M[(i+1):I,1:j] ),
                                           sum( M[1:i,(j+1):J] ),
                                           sum( M[(i+1):I,(j+1):J] ) ),
                                        nrow = 2 ) )
    }
  }
  return(OR_global)
}


local_OR( DoseResult )

global_OR( DoseResult )

ffold_local <- function ( data ){
  # I and J
  I <- nrow(data)
  J <- ncol(data)
  par( mfrow = c(I-1, J-1) )
  for( i in 1:(I-1) ){
    for( j in 1:(J-1) ){
      sub_data <- data[c(i,i+1),c(j,j+1)]
      fourfoldplot( sub_data )
    }
  }
}

# Apply fourfold local OR plot function.
ffold_local( DoseResult )

score_x <- 1:3
score_y <- 1:3
score_xy <- crossprod( t(score_x), t(score_y) )
n <- sum(DoseResult)
sum_xy <- sum( score_xy * DoseResult )
sum_x <- sum( score_x * rowSums(DoseResult) )
sum_y <- sum( score_y * colSums(DoseResult) )
sum_x2 <- sum( score_x^2 * rowSums(DoseResult) )
sum_y2 <- sum( score_y^2 * colSums(DoseResult) )
rxy <- ( n * sum_xy - sum_x * sum_y ) / ( sqrt( n * sum_x2 - sum_x^2 ) * sqrt( n * sum_y2 - sum_y^2 ) )
M2 <- (n-1) * rxy^2
p.value <- 1-pchisq(M2,1)
linear.trend( table = DoseResult, x = score_x, y = score_y )



Titanic
dim( Titanic )
dimnames( Titanic )

Titanic[,, Age="Child", Survived="No"]
margin.table( Titanic, margin = c(1,2) )
margin.table( Titanic, margin = 4 )
nominal_OR( margin.table( Titanic, margin = c(1,2) ) )
local_OR( margin.table( Titanic, margin = c(1,2) ) )
global_OR( margin.table( Titanic, margin = c(1,2) ) )
chisq.test( margin.table( Titanic, margin = c(1,4) ) )

# An unsurprising overwhelming amount of evidence
# that these two are associated.

linear.trend( table = margin.table( Titanic, margin = c(1,4) ),
              x = 1:4, y = 1:2 )
sieve( Titanic )
mosaic( Titanic )
