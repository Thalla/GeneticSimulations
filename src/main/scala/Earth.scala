import AA._

import scala.util.Random
import java.io.{BufferedWriter, File, FileWriter}

import com.github.tototoshi.csv._

import scala.collection.mutable.{ArrayBuffer, ListBuffer}  //TODO?: replace by own methods (data-> List-> needed type => data-> needed type)


/** Fourth draft
  *
  * @nextSteps: similarity, probabilities, NN?, DM?, evo. Alg.?, analyse, parameters, Empirie, Fitness, Metrik
  * @author Hanna Schumacher
  * @version 4.1   -> 500.000 in 315s ; in 2?? s; in 272s without some stuff, 190s without all
  * TODO choose by translation fitness, not always the best one, use pathvariable
  */
object Earth extends Cell{
  //TODO: check if start parameters are valid
  //TODO: tests, exceptions, logging, class diagram, cases, directly connect code/model decisions to research
  /** initiates and starts simulation
    * @param args
    */







  /*
   * simulation; holds global parameters; cares for cell reproduction
   */
  def simulate(cell:Cell, steps:Int):Array[Array[Int]] = {



    /*



    newComparison = nextCell.aaHasTransl
      // add line to aaChanges
      aaChanges(nextCell.generationID) = new Array[Int](2)
      for(
        i <- 0 until aaNumb
      ){
        (oldComparison(i), newComparison(i)) match {
          case (false, true) => {
            aaChanges(nextCell.generationID - 1)(0) += 1
          }
          case (true, false) => aaChanges(nextCell.generationID - 1)(1) += 1
          case _ =>
        }
      }
      oldComparison = newComparison







    file = new File(s"C:\\Users\\feroc\\OneDrive\\Dokumente\\HS\\Semester\\4\\Thesis\\Modeling\\csv\\generationFitness.csv")
    var writer = CSVWriter.open(file)
    writer.writeAll(List(generationFitness))
    writer.close()

    file = new File(s"C:\\Users\\feroc\\OneDrive\\Dokumente\\HS\\Semester\\4\\Thesis\\Modeling\\csv\\mRNAdata.csv")
    writer = CSVWriter.open(file)
    writer.writeAll(mRNAdata)
    writer.close()

    file = new File(s"C:\\Users\\feroc\\OneDrive\\Dokumente\\HS\\Semester\\4\\Thesis\\Modeling\\csv\\testdata.csv")
    bw = new BufferedWriter(new FileWriter(file))
    bw.write(aaChanges.map(_.mkString(", ")).mkString("\n"))
    bw.close()
*/


/*
    //var fitnessCSV:List[List[Double]] = List()

    /*for(
      generation:Array[Double] <- aaNumb
    ){
      fitnessCSV = fitnessCSV :+ generation.toList
    }*/
    file = new File(s"C:\\Users\\feroc\\OneDrive\\Dokumente\\HS\\Semester\\4\\Thesis\\Modeling\\csv\\aaNumb.csv")
    var writer = CSVWriter.open(file)
    writer.writeAll(List(aaNumbPerGen.toList))
    writer.close()







    //aaChanges to file
    var content1 =""
    var content2 =""
                            /*for(
                              i <- 0 until steps
                            ){
                              content1 += i + " "
                              content2 += i + " "
                            }
                            content1 += "\n"
                            content2 += "\n"*/
    /*for(
      i <- 0 until steps
    ){
      content1 += aaChanges(i)(0).toString() + ", "
      //TODO content2 += tableFieldOverTime(i)(0).toString() + ", "
    }
    content1 = content1.dropRight(2)
    content1 += "\n"
    for(
      i <- 0 until steps
    ){
      content1 += aaChanges(i)(1).toString() + ", "
    }
    content1 = content1.dropRight(2)
    content2 = content2.dropRight(2)
    content2 += "\n"
    for(
      i <- 0 until steps
    ){
      //TODO content2 += tableFieldOverTime(i)(1).toString() + ", "
    }
    content2 = content2.dropRight(2)
    content2 += "\n"
    for(
      i <- 0 until steps
    ){
      //TODO content2 += tableFieldOverTime(i)(2).toString() + ", "
    }
    content2 = content2.dropRight(2)
    content2 += "\n"
    for(
      i <- 0 until steps
    ){
      //TODO content2 += tableFieldOverTime(i)(3).toString() + ", "
    }
    content2 = content2.dropRight(2)
    content2 += "\n"
    for(
      i <- 0 until steps
    ){
      //TODO content2 += tableFieldOverTime(i)(4).toString() + ", "
    }
    content2 = content2.dropRight(2)
    val file1 = new File(s"C:\\Users\\feroc\\OneDrive\\Dokumente\\HS\\Semester\\4\\Thesis\\Modeling\\csv\\testdata.csv")
    val bw1 = new BufferedWriter(new FileWriter(file1))
    bw1.write(content1)
    bw1.close()
    val file2 = new File(s"C:\\Users\\feroc\\OneDrive\\Dokumente\\HS\\Semester\\4\\Thesis\\Modeling\\csv\\testdataFirstCell.csv")
    val bw2 = new BufferedWriter(new FileWriter(file2))
    bw2.write(content2)
    bw2.close()*/
*/
    //aaChanges
    Array()
  }




  def writeToFile(filename:String, content:String): Unit ={
    val file = new File(s"C:\\Users\\feroc\\OneDrive\\Dokumente\\HS\\Semester\\4\\Thesis\\$filename.html")
    val bw = new BufferedWriter(new FileWriter(file))
    bw.write(content)
    bw.close()
  }







  /** method that is used for getting random elements with the assurance that no element is chosen twice before all elements have been chosen at least once.
    *
    * @param n a sequence of elements
    * @param nReduced same type as n but having less elements
    * @return a random element either chosen from nReduced or else from n
    */
  def getUniqueRandElem(n:Seq[Any], nReduced:Seq[Any]):Any ={
    val r = new scala.util.Random(22)
    if (nReduced.length >= 1) {
      //get random Tuple
      nReduced(r.nextInt(nReduced.length))
    }
    else { //all elements have been used once, so they can now be used more often
      n(r.nextInt(n.length))
    }
  }

  /** One mRNA holds genes for all initially existing aaRS (see aaRSnumb).
    * Gene length is same as aaRS length. No gene exists twice.
    *
    * @param codons codons/twoTuples from which is randomly chosen to build the mRNA
    * @param geneLength number of codons/twoTuples that code for one aaRS
    * @param geneNumb number of aaRS
    * @return mRNA having genes having codons; for example: Vector(Vector(AA, CU, AG), Vector(UA, CC, GC),...)
    */
  def getRandomMRNA(codons:List[Any], geneLength:Int, geneNumb:Int):List[List[Int]] = {
    var mRNA:List[List[Int]] = List()
    val r = new scala.util.Random(22)
    val getRandomGene = () => {
      var gene:List[Int] = List()
      for {
        _ <- 0 to (geneLength-1)
      } gene = gene :+ r.nextInt(codons.length)
      gene
    }
    // get a gene that isn't already existing
    val getUniqueRandomGene = () => {
      var gene:List[Int] = getRandomGene()
      while(mRNA.contains(gene)){
        gene = getRandomGene()
      }
      gene
    }
    //create mRNA with genes for desired start number of aaRS
    for {
      _ <- 0 to (geneNumb-1)
    } mRNA = mRNA :+ getUniqueRandomGene()

    mRNA
  }

  /**
    * not random
    * @param codons
    * @param geneLength
    * @param geneNumb
    * @param file
    */
  def writeMRNAtoFile(codons:List[Any], geneLength:Int, geneNumb:Int, file:File):Unit = {
    var mRNA:List[List[String]] = List()
    val writer = CSVWriter.open(file)
    val r = new scala.util.Random(22)
    var codonCounter = 0
    for(
      _ <- 0 until geneNumb
    ){
      var gene:List[String] = List()
      for {
        _ <- 0 to (geneLength-1)
      } {
        if(codonCounter < codons.length) {
          gene = gene :+ codonCounter.toString()
          codonCounter += 1
        }else{ //each codon was used once
          gene = gene :+ r.nextInt(codons.length).toString()
        }
      }
      mRNA = mRNA :+ gene
    }

    writer.writeAll(mRNA)
    writer.close()
  }

  /**
    *
    * @param codons
    * @param geneLength
    * @param geneNumb
    * @param file
    */
  def writeAARStoFile(codons:List[Any], aarsLength:Int, aaRSnumb:Int, codonNumb: Int, initAA:Vector[AA], file:File):Unit = {
    val writer = CSVWriter.open(file, append= true)
    val r = new scala.util.Random(22)
    // The current genetic code doesn't have more than six different codons per amino acid. The start cell has the same restrictions.
    val numbOfAnticodonsForOneAARS = List(1,2,3,4,5,6)  //TODO remove redundancy: this is the same as maxAnticodonNumb in Cell

    //create 20^3 aaRS
    for(
      i <- 0 until initAA.length;
      j <- 0 until initAA.length;
      k <- 0 until initAA.length
    ){
      val aaRSline:String = ""
      //var translations:Map[(AA, Int),List[(Double, Int)]] = Map() //each aaRS has a translation table that maps amino acids and codons to a probability and tRNA stems. Currently the probability is not used and always 1.0.
      var antiCodonNumb = numbOfAnticodonsForOneAARS(r.nextInt(numbOfAnticodonsForOneAARS.length))      //number of anticodons that are recognized by the aaRS is chosen randomly

      for(
        _ <- 0 until antiCodonNumb
      ){
        // aa1 aa2 aa3, aa, anticodon, probability, stem
        writer.writeRow(Seq(i, j, k, r.nextInt(initAA.length), r.nextInt(codons.length), r.nextDouble().toString(), r.nextInt(codons.length) ))
      }
    }
    writer.close()
  }

  /**
    * create living aaRS file (create random aaSeqences that will come to life)
    * @return
    */
  def writeLivingAarsToFile():Unit ={
    val writer = CSVWriter.open(new File(path+"livingAARS.csv"), append = true)
    for(
      _ <- 0 until aaRSnumb
    ){
      writer.writeRow(Seq(r.nextInt(aaNumb), r.nextInt(aaNumb), r.nextInt(aaNumb)))
    }
    writer.close()
  }


  def printMRNA(mRNA:List[List[Int]], codons:IndexedSeq[(Any)]):Unit = {
    print(s"mRNA: List of genes. Gene is a List of codon IDs.")
    mRNA.foreach(gene => {
      println()
      gene.foreach(codonID => print(codons(codonID)))
    })
    println()
    println()
  }



  def print2DimTRNAmatrix (matrix:Array[Array[Any]]):Unit = {
    for (
      i <- 0 until matrix.length;
      j <- 0 until matrix(0).length
    ) {
      val elem = matrix(i)(j)
      if (elem == null) print(" ")//i + ", " + j + " ")
      else println(matrix(i)(j).toString())
    }
  }




  def listToString(list:List[Any]):String ={
    var listAsString = ""
    list.foreach(elem => listAsString += elem + " ,")
    listAsString
  }

  // prints cell list
    def toString(cells: Vector[Cell]): String ={
      var toPrint:String = ""
      for(cell <- cells) toPrint += cell.toString()
      toPrint
    }

}