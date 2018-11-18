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
  */
object Earth{
  //TODO: check if start parameters are valid
  //TODO: tests, exceptions, logging, class diagram, cases, directly connect code/model decisions to research
  /** initiates and starts simulation
    * @param args
    */
  def main(args: Array[String]):Unit ={
   /* val steps = 500000
    val aaNumbPerGen = new Array[Int](steps)
    val generationFitness= new Array[Double](steps)
    val aaChanges:Array[Array[Int]] = new Array[Array[Int]](steps)//Array.fill(steps)(Array.fill(2)(0)) //genIds, added amino acids, removed amino acids
    for(
      i <- 0 until 500000
    ){
      aaChanges(i) = new Array[Int](2)
    }
    for(
      i <- 0 until 500000
    ){
      aaChanges(i)(0) = 2
    }
    for(
      i <- 0 until 500000
    ){
      generationFitness(i) = 1.0
    }
    for(
      i <- 0 until 500000
    ){
      aaNumbPerGen(i) = 1
    }
    val tableFieldOverTime:Array[Array[Int]] = new Array[Array[Int]](steps)
    for(
      i <- 0 until 500000
    ){
      tableFieldOverTime(i) = new Array[Int](8)
    }
    for(
      i <- 0 until 500000;
      j <- 0 until 8
    ){
      tableFieldOverTime(i)(j) = 4
    }

    val unambiguousness:Array[Double] = new Array[Double](steps*64)
    for(
      i <- 0 until 500000*64
    ){
      unambiguousness(i) = 5.9
    }
    var x = 5
    var z = 0.0
    for(
      i <- 0 until 500000*64
    ){
      z += (1.0/x.toDouble)/ 63.0 //unambiguousness.foldLeft(0.0)(_+_)
    }

//println(x)*/

   simulate(init(), 500000)
  }


def getX():Array[Array[Int]]={

  val x:Array[Array[Int]] = simulate(init(),1000)
  x
}

  /** creates a start cell, each cell is one generation
    * @param aaRSnumb Number of how many different aaRS proteins exist initially. (see old scaladoc for notes)
    * @param aarsLength The number of amino acids an aaRS consists of. TODO: aaRS can have varying lengths
    * @param initAA initially existing and used amino acids TODO: allow having new/non proteinogenic amino acids TODO start with non essential amino acids
    * @param codonNumb Either 16 (set of twoTuples is created) or 64 (set of real codons is created) or 48 (set of strong commafree codons is created) (this version uses 64, others are not tested)
    * @return start cell
  */
  def init (aaRSnumb:Int = 22, aarsLength:Int = 3, initAA:Vector[AA] = AA.values.toVector, codonNumb:Int = 64):Cell= {
    // create codons
    val codons = getCodons(codonNumb)
    val path = "C:\\Users\\feroc\\OneDrive\\Dokumente\\HS\\Semester\\4\\Thesis\\Modeling\\csv\\"
    val aaNumb = initAA.length

    // two mRNA creation versions, first random with no gene twice, second saved in file
     // create mRNA (List of genes. Each gene is a List of codon IDs as long as aarsLength)
     // val mRNA:List[List[Int]] = getRandomMRNA(codons.toList, aarsLength, aaRSnumb)
     // writeMRNAtoFile(codons.toList, aarsLength, aaRSnumb, new File(path+"mRNA.csv"))

    // create aaRS file
     // writeAARStoFile(codons.toList, aarsLength, aaRSnumb, codonNumb, initAA, new File(path+"aaRS.csv"))

    // create living aaRS file (create random aaSeqences that will come to life)
    /*
    val writer = CSVWriter.open(new File(path+"livingAARS.csv"), append = true)
    for(
      _ <- 0 until aaRSnumb
    ){
      writer.writeRow(Seq(Random.nextInt(aaNumb), Random.nextInt(aaNumb), Random.nextInt(aaNumb)))
    }
    writer.close()
    */


    // read mRNA
    var reader = CSVReader.open(new File(path+"mRNA.csv"))
    var data = reader.all()
    var mRNA:List[List[Int]] = List()
    for(
      line <- data
    ){
      var gene:List[Int] = List()
      for(
        elem <- line
      ){
        gene = gene :+ elem.toInt
      }
      mRNA = mRNA :+ gene
    }
    reader.close()


    // read aaRS
    reader = CSVReader.open(new File(path+"aaRS.csv"))
    data = reader.all()
    var aaSeq:Vector[Int] = Vector()      //position in allAARS
    var translations:Map[(AA, Int),List[(Double, Int)]] = Map()
    // Array with all existing aaRS, initialised with placeholders
    var allAARS:Array[Array[Array[AARS]]] = Array.fill[Array[Array[AARS]]](aaNumb)(Array.fill[Array[AARS]](aaNumb)(Array.fill[AARS](aaNumb)(new AARS(Vector(), translations))))
    // give each aaRS in allAARS the correct aaSeq (is dependent from its position in allAARS)
    for(
      i <- 0 until aaNumb;
      j <- 0 until aaNumb;
      k <- 0 until aaNumb
    ){
      allAARS(i)(j)(k).aaSeq = Vector(initAA(i), initAA(j), initAA(k))
    }
    // read file data and give each aaRS its translations
    for(
      line <- data
    ){
      //get Array Index of aaRS
      val i1 = line(0).toInt
      val i2 = line(1).toInt
      val i3 = line(2).toInt
      // add translation to aaRS
      allAARS(i1)(i2)(i3).translations += ((initAA(line(3).toInt), line(4).toInt)->List((line(5).toDouble, line(6).toInt)))
    }
    reader.close()


    // read the living aaRS
    var livingAARSs:List[AARS] = List()
    val lifeticksStartValue = allAARS(0)(0)(0).lifeticksStartValue
    reader = CSVReader.open(new File(path+"livingAARS.csv"))
    data = reader.all()
    for(
      line <- data
    ){
      allAARS(line(0).toInt)(line(1).toInt)(line(2).toInt).lifeticks = lifeticksStartValue
      livingAARSs = livingAARSs :+ allAARS(line(0).toInt)(line(1).toInt)(line(2).toInt)
    }
    reader.close()


    //create code Table: Array[Codon][AA] = List[AARS]
    val codeTable:Array[Array[List[AARS]]] = Array.ofDim[List[AARS]](codons.length, aaNumb)
    livingAARSs.foreach(aaRS => {
      aaRS.translations.foreach(translation => {
        var pos = codeTable(translation._1._2)(translation._1._1.id)
          if(pos == null){
            codeTable(translation._1._2)(translation._1._1.id) = List(aaRS)
          } else{
            codeTable(translation._1._2)(translation._1._1.id) = codeTable(translation._1._2)(translation._1._1.id) :+ aaRS
          }
      })
    })

    new Cell(mRNA, codeTable, livingAARSs.toVector, allAARS, initAA, codonNumb,0)

  }




  /*
   * simulation; holds global parameters; cares for cell reproduction
   */
  def simulate(cell:Cell, steps:Int):Array[Array[Int]] = {//: Vector[Cell] = {
    //var buildUpContent = new StringBuilder
    //buildUpContent ++= cell.toHtmlString(List(PrintElem.mRNA, PrintElem.livingAARSs, PrintElem.codeTable))
    /*val newComparison = new Array[Boolean](cell.initAA.length)  // true if aa has translation
    var oldComparison = new Array[Boolean](cell.initAA.length) //Array.fill(20)(false)
    var numbTranslatedAA = 19 //TODO
    val aaNumbPerGen = new Array[Int](steps)
    val generationFitness= new Array[Double](steps)      //Array.fill(steps)(Array.fill(2)(0.0))
    val aaNumb = cell.initAA.length
    val codonNumb = cell.codonNumb
    //TODOvar mRNAdata:Array[List[String]] = new Array[List[String]](steps*codonNumb)//Array.empty[List[String]]


    //initiate oldComparison
    for( //foreach aa
      i <- 0 until aaNumb
    ){ var goOn = true //stop as soon as a translation for the aa is found
      for(
        j <- 0 until codonNumb if goOn != false
      ){
         if(cell.codeTable(j)(i) != null){
          goOn = false
          oldComparison(i) = true
        }

      }
    if(goOn == true) oldComparison(i) = false
    }

    //TODO val tableFieldOverTime:Array[Array[Int]] = new Array[Array[Int]](steps)//Array.fill(steps)(Array.fill(8)(0))  //assumption: max numb of aaRS per field is 8 TODO
    var value = 0
    if(cell.codeTable(0) != null && cell.codeTable(0)(0) != null){
      value = (cell.codeTable(0)(0)(0)).aaSeq(0).id + (cell.codeTable(0)(0)(0)).aaSeq(1).id*20 + (cell.codeTable(0)(0)(0)).aaSeq(2).id*40 //value is calculated aaRS id
    }
    //TODO tableFieldOverTime(cell.generationID) = new Array[Int](8) //TODO
    //TODO tableFieldOverTime(cell.generationID)(0)= value

    val aaChanges:Array[Array[Int]] = new Array[Array[Int]](steps)//Array.fill(steps)(Array.fill(2)(0)) //genIds, added amino acids, removed amino acids
    aaChanges(0) = new Array[Int](2)
    */
    var newCell = cell
    def time[R](block: => R): R = {
      val t0 = System.nanoTime()
      val result = block    // call-by-name
      val t1 = System.nanoTime()
      println("Elapsed time: " + (t1 - t0) + "ns")
      result
    }

    time(while(newCell.generationID <= steps-2){

      val x = newCell.translate()
      /*    //fitness(newCell.generationID)(1) = newCell.unambiguousness/66
      aaNumbPerGen(newCell.generationID)= numbTranslatedAA
          //newCell.unambiguousness = newCell.unambiguousness.dropRight(2)
      generationFitness(newCell.generationID)= (numbTranslatedAA.toDouble/aaNumb.toDouble)* newCell.unambiguousness   //((newCell.unambiguousness.foldLeft(0.0)(_+_)) /codonNumb.toDouble
      numbTranslatedAA = 0
      for(
        i <- 0 until codonNumb
      ){
        //TODO mRNAdata(newCell.generationID*codonNumb+codonNumb+i) = List(newCell.generationID.toString(), i.toString(), newCell.mRNAdata(i))
      }
      //val row:List[String] = newCell.mRNAdata
      //TODO mRNAdata = mRNAdata :+ row*/
      newCell = x

      if(newCell.generationID%100000 == 0){
        val x = "hmm"//buildUpContent ++= newCell.toHtmlString(List(PrintElem.mRNA, PrintElem.livingAARSs, PrintElem.codeTable))
      }

    /*  //TODO tableFieldOverTime(newCell.generationID) = new Array[Int](8) //TODO
      if(newCell.codeTable(0)(0) != null){
        val z = newCell.codeTable(0)(0).length
        for(
          i <- 0 until newCell.codeTable(0)(0).length
        ){
          if(newCell.codeTable(0)(0)(i) != null){
            value = (newCell.codeTable(0)(0)(i)).aaSeq(0).id + (newCell.codeTable(0)(0)(i)).aaSeq(1).id*20 + (newCell.codeTable(0)(0)(i)).aaSeq(2).id*40
            //TODO tableFieldOverTime(newCell.generationID)(i) = value
          }
        }
      }


      //update newComparison
      for(
        i <- 0 until aaNumb  // go through amino acids TODO
      ){
          var goOn = true
          for(
          j <- 0 until codonNumb if goOn != false  //TODO
        ){
          if(newCell.codeTable(j)(i) != null){
            goOn = false
            newComparison(i) = true
          }
        }
        if(goOn == true) newComparison(i) = false
      }
      // add line to aaChanges
      aaChanges(newCell.generationID) = new Array[Int](2)
      for(
        i <- 0 until aaNumb
      ){
        (oldComparison(i), newComparison(i)) match {
          case (false, true) => {
            aaChanges(newCell.generationID - 1)(0) += 1
            numbTranslatedAA += 1
          }
          case (true, false) => aaChanges(newCell.generationID - 1)(1) += 1
          case (true, true) => numbTranslatedAA += 1
          case _ =>
        }
      }
      oldComparison = newComparison
      */
      //println(newCell.generationID)
    }

    )
   /* val content = "HMM"//buildUpContent + newCell.toHtmlString(List(PrintElem.mRNA, PrintElem.livingAARSs, PrintElem.codeTable))
    var file = new File(s"C:\\Users\\feroc\\OneDrive\\Dokumente\\HS\\Semester\\4\\Thesis\\Modeling\\Scala\\simulationOutput.html")
    val bw = new BufferedWriter(new FileWriter(file))
    bw.write(content)
    bw.close()

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


    file = new File(s"C:\\Users\\feroc\\OneDrive\\Dokumente\\HS\\Semester\\4\\Thesis\\Modeling\\csv\\generationFitness.csv")
    writer = CSVWriter.open(file)
    writer.writeAll(List(generationFitness.toList))
    writer.close()

    file = new File(s"C:\\Users\\feroc\\OneDrive\\Dokumente\\HS\\Semester\\4\\Thesis\\Modeling\\csv\\mRNAdata.csv")
    writer = CSVWriter.open(file)
    //writer.writeAll(mRNAdata)
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


  /** generate from nucleobases
    * @param codonNumb 16 -> twoTuples, 64 -> codons, 48 -> strong comma free codons
    * @return either list of twoTuples or codons or strong comma free codons
    */
  def getCodons(codonNumb:Int):IndexedSeq[(Any)] = {
    //range for adding third base
    val range:Range = codonNumb  match{
      case 16 => 0 to 0       //-> no third base
      case 64 | 48 => 0 to 3
    }
    // if strong comma free code shall be generated this filter is used
    val filterScf = (k:Int, i:Int) => {
      if(codonNumb == 48) k!=i
      else true
    }

    var codons = (for {
      i <- 0 to 3 // -> four nucleobases
      j <- 0 to 3
      k <- range
      if filterScf(k, i)  // if scf option 48 then filter codons
    } yield {
      codonNumb match{
        case 16 => (Base(i),Base(j))
        case 64 | 48=> (Base(i),Base(j), Base(k))
      }})
    codons
  }




  /** method that is used for getting random elements with the assurance that no element is chosen twice before all elements have been chosen at least once.
    *
    * @param n a sequence of elements
    * @param nReduced same type as n but having less elements
    * @return a random element either chosen from nReduced or else from n
    */
  def getUniqueRandElem(n:Seq[Any], nReduced:Seq[Any]):Any ={
    if (nReduced.length >= 1) {
      //get random Tuple
      nReduced(Random.nextInt(nReduced.length))
    }
    else { //all elements have been used once, so they can now be used more often
      n(Random.nextInt(n.length))
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

    val getRandomGene = () => {
      var gene:List[Int] = List()
      for {
        _ <- 0 to (geneLength-1)
      } gene = gene :+ Random.nextInt(codons.length)
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
          gene = gene :+ Random.nextInt(codons.length).toString()
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

    // The current genetic code doesn't have more than six different codons per amino acid. The start cell has the same restrictions.
    val numbOfAnticodonsForOneAARS = List(1,2,3,4,5,6)

    //create 20^3 aaRS
    for(
      i <- 0 until initAA.length;
      j <- 0 until initAA.length;
      k <- 0 until initAA.length
    ){
      val aaRSline:String = ""
      //var translations:Map[(AA, Int),List[(Double, Int)]] = Map() //each aaRS has a translation table that maps amino acids and codons to a probability and tRNA stems. Currently the probability is not used and always 1.0.
      var antiCodonNumb = numbOfAnticodonsForOneAARS(Random.nextInt(numbOfAnticodonsForOneAARS.length))      //number of anticodons that are recognized by the aaRS is chosen randomly

      for(
        _ <- 0 until antiCodonNumb
      ){
        // aa1 aa2 aa3, aa, anticodon, probability, stem
        writer.writeRow(Seq(i, j, k, Random.nextInt(initAA.length), Random.nextInt(codons.length), Random.nextDouble().toString(), Random.nextInt(codons.length) ))
      }
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

  def mRNAtoHTML(mRNA:List[List[Int]], codons:IndexedSeq[(Any)]):String = {
    var content = "<ul>\n"
    mRNA.foreach(gene => {
      content += "<li>\n"
      gene.foreach(codonID => content += codons(codonID).toString() )
      content += "</li>\n"
    })
    content += "</ul>\n"
   content
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

  //print aaRS-Matrix
  def print3DimAARSmatrix (matrix:Array[Array[Array[AARS]]]):Unit = {
    for (
      i <- 0 until matrix.length;
      j <- 0 until matrix(0).length;
      k <- 0 until matrix(0)(0).length
    ) {
      val elem = matrix(i)(j)(k)
      if (elem == null) {}//println(i + ", " + j + ", " + k + " ")
      else println(i + " , " + j + " , " + k + ": " + matrix(i)(j)(k).toString())
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