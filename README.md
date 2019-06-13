# GeneticSimulations
The code for analysing the simulation results can be found in repository [GeneticSimulationsR](https://github.com/Thalla/GeneticSimulationsR).  

In the case of questions and/or comments, which are highly appreciated, write to temerita16@gmail.com .

Running the main method of object Simulator with an empty args parameter results in a couple of simulations and the results are saved in the subfolder "SimulationResults". 

**Caution:** When simulation breaks, all folders of files that may not be finished correctly must be deleted before running the simulation again. 
When a file shall be recreated, the folder of the file must be deleted.


## Parameter descriptions
- **basePath:** path were the folders with the results of the simulation are created
- **codonNumb:** number of codons that is used for mRNA and aaRSs translations creation
- **livingAarsSeed:** seed for random number generator that decides on which living aaRSs the simulation starts with
- **steps:** number of generations that shall be simulated
- **translationMethod:** defines how the decision which translation will be taken is made in the case of an unambiguous codon  
  </br>
  
- **mrnaSeed:** seed value for the random number generator that chooses the codons used for the mRNA after all codons have been used once
- **geneLength:** defines length of aaRS sequence and is always three
- **geneNumbs:** number of genes that are translated in each cell generation
- **similarAarss:** true -> similar aaRSs sequences also have similar translations
- **initAaNumbs:** number of amino acids that are available in the simulation, maximal 23, for more amino acids Enum AA must be extended
- **maxAnticodonNumbs:** maximal number of translations an aaRSs can have
- **aarsLifeticksStartValues:** defines how long aaRSs live. An aaRS with ten lifeticks survives ten generations.
- **aarsSeed seed value:** for random number generator that defines the translations each aaRS has
- **livingAarsStartNumbs:** number of living aaRSs in the beginning of the simulation
- **outputSeed:** seed for random number generator that influences which translation is taken in the case of translation method "random" or "affinity"
- **addToOutput:** true -> the last set of living aaRSs is saved in a file. A second run of the same parameter combinations starts with the last saved set of living aaRSs.


## Main Objects

### Simulator
**Parameters for main method:**
- basePath
- codonNumb
- livingAarsSeed
- steps
- translationMethod

**Predefined parameters:**
- mrnaSeed: 2000
- geneLength: 3
- geneNumbs: 2, 5, 15, 20, 30
- similarAarss: false, true
- initAaNumbs: 3, 4, 8, 10, 16, 20
- maxAnticodonNumb: 2, 4, 6, 8
- aarsLifeticksStartValue: 2, 5, 10, 15
- aarsSeed: 200
- livingAarsStartNumbs: 2, 5, 15, 20, 30
- outputSeed: 1
- addToOutput: true


These are the default values of the initSimulation() method.


### Reduced Simulator
**Parameters for main method:**
- basePath
- codonNumb
- livingAarsSeed
- steps
- translationMethod
- initAaNumbs
- outputSeed

**Predefined parameters**:
- geneNumbs: 5, 10, 20, 30
- similarAarss: false, true
- maxAnticodonNumb: 6
- aarsLifeticksStartValue: 2, 5, 10, 15
- livingAarsStartNumbs: 5, 10, 20, 30


The remaining parameters are used with their default values (see Simulator).


### Mini Simulator
**Parameters for main method:**
- basePath
- codonNumb
- livingAarsSeed
- steps
- translationMethod

**Predefined parameters:**
- geneNumbs: 1, 2, 3, 4, 5
- similarAarss: false, true
- initAaNumbs: 1, 2, 3, 4, 5
- maxAnticodonNumbs: 1, 2, 3, 4, 5, 6, 7, 8
- aarsLifeticksStartValues: 1, 2, 3, 4, 5, 6, 7, 8, 9, 10
- livingAarsStartNumbs: 1, 2, 3, 4, 5


The remaining parameters are used with their default values (see Simulator).


### Custom Simulator
The main method needs values for all parameters.

# Files
**Naming Convention:**  
*SimulationResults\mRNA_s200_c64_gl3_gn20\initAaNumb_20\similar_false\aaRS_s30_ac6_lt20\livingAars_n10_s20\output_s1*
This path leads to the fitness file of a simulation with the following parameters:
- mRNA with seed value 200, codon number 64, gene length 3, gene number 20
- 20 amino acids
- aaRSs with no similarity
- aaRSs with seed value 30, anticodon number 6 and 20 lifeticks
- 10 living aaRSs with seed value 20
- output seed 1

# How to make changes
**Change the simulation method**  
Currently the simulation method is *simulateCell()*. *runSimulation()* needs a simulation method as input as well as a list of cells. The cells are simulated in parallel with the simulation method. *simulateCell* triggers a simulation step and afterwards the update of *SimulationData*, which writes data as the fitness values to the files, till the desired number of steps is reached. *simulateCell* can be replaced with a method that handles the simulation differently, for example:
- Initialize more than one cell. For example: cell1 = Cell(outputSeed); cell2 = Cell(outputSeed+1)
- Simulate natural selection: Take only the cells with the hightest fitness values and with a given probability cells with a lower fitness value.
- Simulate cell/membrane division: Divide the set of living aaRSs. This can lead to much higher fitness values because less aaRSs and therefor less translations allows for less ambiguity.
- The fittest wins: A translation step is simulated with various seeds and the fittest result is the cell for the next generation.

**Simulate special cases**  
To simulate a specific set of genes or aaRSs or living aaRSs the simulation must be run with the correct parameters first and afterwards the desired detail changes can be made in the files that were generated. A second run of the simulation will use the changed files and perform the desired simulation. A second run is only possible when the ouput folder was deleted in before. Otherwise, when the first run had addToOutput = false nothing happens and when addToOutput was true the second simulation uses the last set of aaRSs from the first simulation.

# Future Features
**Feedback on which feature shall be implemented next is highly appreciated**
- Cell division: set of aaRS is divided into two sets. This leads to an overall fitness improvement.
  - Two or more possible divisions are made, the case that leads to the best fitness is used for furhter simulation.
- Not all peptides with length three are valid aaRSs. 
- Peptides that are no aaRSs have some kind of functionality and/or a positive or negative influence on translation fidelity etc.
- mRNA mutations
  - random
  - directed so that fitness improves
