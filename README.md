# Tools
A repository to hold generalizable functions and tools for use by lab members and collaborators

## Usage
To use any of the Crombie Lab R functions you can source the `code/functions.R` file directly from this repository from RStudio. Just run this line of code
```
eval(parse(text = readLines("https://github.com/Crombie-Lab/Tools/raw/main/code/functions.R")))
```

## Function help
---
### geneticDistance()
This function will calculate the genetic distance between two physical positions in the C. elegans genome. The arguments to supply are `left` and `right` physical positions and `chrom`. Here's an example of how to use the function and its output.
```
g.dist <- geneticDistance(left = 3184137, right = 4436190, chrom = "V")
getting C. elegans genetic map from https://github.com/AndersenLab/post-gatk-nf/raw/main/input_files/annotations/c_elegans/c_elegans_genetic_map.bed.gz
 [100%] Downloaded 170878 bytes...
The genetic distance between V:3184137 and V:4436190 is 7.3644
```
