# Practical 1 - Contingency tables


DR_data <- matrix( c(41, 9,
                     37, 13 ), byrow = TRUE, ncol = 2 )

dimnames( DR_data ) <- list( Dose = c("High", "Low"),
                             Result = c("Success", "Failure") )

DR_data


DR_contingency_table <- addmargins( DR_data )
DR_contingency_table

DR_prop <- prop.table( DR_data )
DR_prop_table <- addmargins( DR_prop )  # adds marginals (sums) on the bottom and right of the table
DR_prop_table


prop.table(DR_data, 1)  # row conditionals


DR_prop_1 <- prop.table( DR_data, 1 )
DR_prop_1_table <- addmargins( DR_prop_1, margin = 2 )
DR_prop_1_table


DR_prop_2_table <- addmargins( DR_prop_1_table, margin = 1, FUN = mean ) # does mean function instead of sum
DR_prop_2_table


DR_prop_3_table <- addmargins( DR_prop_1,
                               margin = c(1,2),
                               FUN = list(mean, sum) )

DR_prop_3_table ## same thing


####### Construction from dataframes


library(palmerpenguins)
head(penguins)
?penguins


penguins_data <- table( Species = penguins$species, Island = penguins$island )
penguins_data



## Exercises

# a
penguins_margins <- addmargins(penguins_data)
penguins_margins

# b
penguins_proportions <- prop.table(penguins_data)
penguins_proportions

# c
penguins_col_cond <-  prop.table(penguins_data, 2)
penguins_col_cond

penguins_col_cond_marg <- addmargins(penguins_col_cond, margin = 1, FUN = sum)
penguins_col_cond_marg

# d

penguin_species_prop <- addmargins(penguins_proportions, margin = 2, FUN = sum)
penguin_species_prop
penguins_final <- cbind(penguins_col_cond_marg, c(penguin_species_prop[,4], 1))
penguins_final




### Chi square test

chisq.test( DR_data, correct = FALSE ) # no continuity correction

chisq.test( DR_data ) # chi^2 test with continuity correction

chisq.test(DR_data)$expected # ML estimates of the expected frequencies


chisq.test(penguins$species, penguins$island)
# p-value = 2.2x10^-16, so it is v significant


### Data Visualisation

DR_prop

barplot(DR_prop)

?barplot

barplot(DR_prop, density = 70)
barplot(DR_prop, density = 30)
barplot(DR_prop, density = 0)

title("Barplot",
      xlab = "success/failure",
      ylab = "proportion")

?barplot

barplot(DR_prop,
        density = 70,
        xlab = "success/failure",
        ylab = "proportion",
        legend.text = TRUE)

title("Barplot")


barplot(prop.table(DR_data, 2),
        density = 70,
        xlab = "result",
        ylab = "proportion",
        legend.text = TRUE)

barplot(t(prop.table(DR_data)),
        density = prop.table(DR_data, 1)[,1]*100,
        xlab = "dose level",
        ylab = "proportion",
        legend.text = TRUE)


fourfoldplot( DR_data )
fourfoldplot( DR_data, color = c("red", "blue") )


library(vcd)
sieve( DR_data )

library(vcd)
sieve( DR_data, shade = T )

sieve( DR_data, sievetype = "expected", shade = T )


mosaic( DR_data )


### Odds Ratios in R
# a

odd <- function(p) {
  odd <- p/(1 - p)
  return(odd)
}

odds_ratio <- function(pA, pB) {
  ratio_AB <- odd(pA)/odd(pB)
  return(ratio_AB)
}


odds_ratio_table <- function(data) {
  r <- data[1,1] * data[2,2] / (data[1,2] * data[2,1])
  return(r)
}

odds_ratio_table(DR_data)

# yes if p_12, or p_21 = 0 it would divide by 0

# yes because then one of htem would be 0

odds_ratio_table_adjusted <- function(data) {
  if (data[1,1] == 0) {
    print("Error message: 1,1 = 0")
    stop()
  } else if (data[1,2] == 0) {
    print("Error message: 1,2 = 0")
    stop()
  } else if (data[2,1] == 0) {
    print("Error message: 2,1 = 0")
    stop()
  } else if (data[2,2] == 0) {
    print("Error message: 2,2 = 0")
    stop()
  } else {
    r <- data[1,1] * data[2,2] / (data[1,2] * data[2,1])
    return(r)
  }
}


odds_ratio_table_adjusted_2 <- function(data) {
  if (data[1,1] == 0) {
    print("Error message: 1,1 = 0")
  } else if (data[1,2] == 0) {
    print("Error message: 1,2 = 0")
  } else if (data[2,1] == 0) {
    print("Error message: 2,1 = 0")
  } else if (data[2,2] == 0) {
    print("Error message: 2,2 = 0")
  } else {
    r <- data[1,1] * data[2,2] / (data[1,2] * data[2,1])
    return(r)
    stop()
  }
  warning("Add 0.5")
  data <- data + 0.5
  r <- data[1,1] * data[2,2] / (data[1,2] * data[2,1])
  return(r)
}


###### 5.1.5 Mushrooms

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
