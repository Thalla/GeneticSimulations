import org.ddahl.rscala._
import java.io.{BufferedWriter, File, FileWriter}

import EncodedAaTracker.Update

import scala.concurrent.duration._
import akka.pattern.ask
import PrintElem.PrintElem
import akka.actor.{ActorRef, ActorSystem, PoisonPill}
import akka.util.Timeout

import scala.concurrent.Await
import scala.swing._
import scala.swing.event.ButtonClicked
//import org.biojava3.core.sequence.DNASequence
import AA._
import Base.Base
import scala.util.Random


/** Third draft
  * tRNAs have a less important role
  * Uses 64 codons and 20 amino acids.
  * The code table holds multiple translations.
  * @nextSteps: similarity, probabilities, NN?, DM?, evo. Alg.?, analyse, parameters, Empirie
  * @author Hanna Schumacher
  * @version 0.6 -> + code table - circle relationship
  */
object Earth{
  //TODO: check if start parameters are valid
  //TODO: tests, exceptions, logging, class diagram, cases, directly connect code/model decisions to research

  /** initiates and starts simulation
    * @param args
    */
  def main(args: Array[String]): Unit ={
    val data = new Data()
    val z = data.randomTable()
  simulate(init())
    //val auth = new InputParamsDialog().auth.getOrElse(throw new IllegalStateException("You should login!!!"))
    //println(auth.toString())
    /*val cells:List[Cell] = simulate(init()).toList
    val data = new Data(cells.toArray)
    val codeTableCountAARS:Array[Array[Int]] = Array.ofDim[Int](cells(0).codonNumb, cells(0).initAA.length)

   val hasTranslationOverTime= for(
      c <- 0 until cells.length
    )yield{
     val codeTableHasTranslation:Array[Array[Int]] = Array.ofDim[Int](cells(0).codonNumb, cells(0).initAA.length)
     for(
       i <- 0 until 64;
       j <- 0 until 20
     ){
       val elem = cells(c).codeTable(i)(j)
       /*val elem2 = cells(c+1).codeTable(i)(j)
       def getLength (e:List[AARS]):Int = {
         if(e == null){
           0
         }else{
           e.length
         }
       }*/
       if(elem != null){
         codeTableHasTranslation(i)(j) = 1
       }else{
         codeTableHasTranslation(i)(j) = 0
       }
       codeTableHasTranslation
     }

codeTableHasTranslation
    }

    val file = new File(s"C:\\Users\\feroc\\OneDrive\\Dokumente\\Thesis\\booleanData.txt")
    val bw = new BufferedWriter(new FileWriter(file))
    for(
      i <- 0 until 64;
      j <- 0 until 20
    ){
      bw.write("\n")
      for(
        k <- 0 until cells.length
      ){
        bw.write(hasTranslationOverTime(k)(i)(j) + ", ")
      }

    }

    bw.close()
    var x = 6+6*/
  }


  /** creates a start cell, each cell is one generation
    * @param aaRSnumb Number of how many different aaRS proteins exist initially. (see old scaladoc for notes)
    * @param aarsLength The number of amino acids an aaRS consists of. TODO: aaRS can have varying lengths
    * @param initAA initially existing and used amino acids TODO: allow having new/non proteinogenic amino acids TODO start with non essential amino acids
    * @param codonNumb Either 16 (set of twoTuples is created) or 64 (set of real codons is created) or 48 (set of strong commafree codons is created) (this version uses 64, others are not tested)
    * @return start cell
  */
  def init (aaRSnumb:Int = 20, aarsLength:Int = 3, initAA:Vector[AA] = AA.values.toVector, codonNumb:Int = 64):Cell= {
    // create codons
    val codons = getCodons(codonNumb)

    // start connection is: codon with id x has tRNA with id x has aaRS with id x (three tRNA have a second aaRS) has amino acid with id x

    // create random mRNA (List of genes. Each gene is a List of codon IDs as long as aarsLength)
    val mRNA:List[List[Int]] = getRandomMRNA(codons.toList, aarsLength, aaRSnumb)

    // prepare creation of aaRS
    val allAARS =  Array.ofDim[AARS](initAA.length,initAA.length, initAA.length)    // Array that saves living and dead aaRS. Their position is defined by their amino acid sequence        TODO ? more dimensions ?
    val aaRSsequences:IndexedSeq[IndexedSeq[AA]] = for(          // prepare random amino acid sequences to provide them later to the aaRS constructor
      _ <- 0 until aaRSnumb
    )yield{
      val seq:IndexedSeq[AA] = for(
        j <- 0 until aarsLength
      ) yield AA(Random.nextInt(initAA.length))
      seq
    }  //TODO: check if any seqence is twice (find method in older code)

    // The current genetic code doesn't have more than six different codons per amino acid. The start cell has the same restrictions.
    val numbOfAnticodonsForOneAARS = List(1,2,3,4,5,6)
    var numbUsedCodons = 0

    //create the start set of aaRS, codon with id 63 isn't used and therefore behaves as stop codon
    val livingAARSs = for(
      i <- 0 until aaRSnumb
    )yield{
      var translations:Map[(AA, Int),List[(Double, Int)]] = Map() //each aaRS has a translation table that maps amino acids and codons to a probability and tRNA stems. Currently the probability is not used and always 1.0.
      var antiCodonNumb = numbOfAnticodonsForOneAARS(Random.nextInt(numbOfAnticodonsForOneAARS.length))      //number of anticodons that are recognized by the aaRS is chosen randomly
      if((aaRSnumb-i) > (codonNumb-numbUsedCodons-antiCodonNumb)){
        antiCodonNumb -= (aaRSnumb-i) - (codonNumb-numbUsedCodons-antiCodonNumb)
      }
      for(
        j <- 0 until antiCodonNumb
      ){
          translations += (initAA(i), numbUsedCodons)->List((1.0, Random.nextInt(codons.length)))
          numbUsedCodons += 1
      }

      new AARS(aaRSsequences(i).toVector, translations)
    }

    //create code Table: Array[Codon][AA] = List[AARS]
    val codeTable:Array[Array[List[AARS]]] = Array.ofDim[List[AARS]](codons.length, initAA.length)
    livingAARSs.foreach(aaRS => {
      aaRS.translations.foreach(translation => {
        codeTable(translation._1._2)(translation._1._1.id) = List(aaRS)
        allAARS(aaRS.aaSeq(0).id)(aaRS.aaSeq(1).id)(aaRS.aaSeq(2).id) = aaRS
      })
    })

    new Cell(mRNA, codeTable, livingAARSs.toVector, allAARS, initAA, codonNumb,0)

  }




  // simulation; holds global parameters; cares for cell reproduction
  def simulate(cell:Cell):Unit = {//: Vector[Cell] = {
    //Each cell represents one generation and Earth holds all.
    //var cells: Vector[Cell] = Vector()
    //cells = cells :+ cell



   /* implicit val timeout: Timeout = Timeout(500 seconds)
    val _system: ActorSystem = ActorSystem.create("ColCal")
    val actor: ActorRef = _system.actorOf(EncodedAaTracker.props(Array()))
    actor ! EncodedAaTracker.Update(1,true)
    actor ! 1
    val x= (actor ? Array()).mapTo[Array[Array[Int]]]
    val xx= Await.result(x, 500 seconds)
    _system.terminate()
*/

    var comparisons:Array[Array[Boolean]] = Array.ofDim[Boolean](3,20) //TODO
    var encodedAaChanges:Array[Array[Int]] = Array.ofDim[Int](1000000,20) //TODO
    comparisons(0) = Array.fill(20){true}
    comparisons(0)(19) = false
    var oldGen = 0
    var newGen = 1

    var newCell = cell
    def time[R](block: => R): R = {
      val t0 = System.nanoTime()
      val result = block    // call-by-name
      val t1 = System.nanoTime()
      println("Elapsed time: " + (t1 - t0) + "ns")
      result
    }
    time(while(newCell.generationID < 500000){
      /*for(
        i <- 0 until 19
      ) {
        val gettf = ()=> {
          for (
            j <- 0 until 64
          ) {
            if (newCell.codeTable(j)(i)(0).isInstanceOf[AARS]) {
              return true
            }
          }
          return false
        }
        comparisons(oldGen)(i) = gettf()

      }*/

      newCell = newCell.translate()
/*
      for(
        i <- 0 until 19
      ) {
        val gettf = ()=> {
          for (
            j <- 0 until 64
          ) {
            if (newCell.codeTable(j)(i)(0).isInstanceOf[AARS]) {
              return true
            }
          }
          return false
        }
        comparisons(newGen)(i) = gettf()

      }


      for(
        i <- 0 until 19
      ){
        (comparisons(oldGen)(i), comparisons(newGen)(i)) match {
          case (false, false) => encodedAaChanges(newCell.generationID - 1)(i) = 0
          case (false, true) => encodedAaChanges(newCell.generationID - 1)(i) = 1
          case (true, false) => encodedAaChanges(newCell.generationID - 1)(i) = -1
          case (true, true) => encodedAaChanges(newCell.generationID - 1)(i) = 0
        }
      }
*/

      //println(newCell.generationID)
    })
    val content = newCell.toHtmlString(List(PrintElem.mRNA, PrintElem.livingAARSs, PrintElem.codeTable))
    val file = new File(s"C:\\Users\\feroc\\OneDrive\\Dokumente\\Thesis\\simulationOutput.html")
    val bw = new BufferedWriter(new FileWriter(file))
    bw.write(content)
    bw.close()




    //finally bw.close()


    //print cell
    //print(cell.toString())
    /*var content = cell.toHtmlString(List(PrintElem.codons, PrintElem.mRNA, PrintElem.livingAARSs, PrintElem.codeTable))
    //writeToFile("generation0", content)
    val file = new File(s"C:\\Users\\feroc\\OneDrive\\Dokumente\\Thesis\\simulationOutput.html")
    val bw = new BufferedWriter(new FileWriter(file))
    bw.write(content)*/


    /*def go(cell:Cell):Unit = {*/

      //val newCell:Cell = cell.translate()
      //cells = cells :+ newCell
      //writeToFile("generation"+cell.generationID, content)
      /*if(!(newCell.livingAARSs.length != 0)) {
        content = newCell.toHtmlString(List(PrintElem.livingAARSs, PrintElem.codeTable))
        bw.write(content)
      }else{
        if(newCell.generationID%100000 == 0){
          //print("GENERATION" + cell.generationID + cell.toString())
          /*content = newCell.toHtmlString(List(PrintElem.livingAARSs, PrintElem.codeTable))
          bw.write(content)*/
        }*/
        //go(newCell)

      //}
   /* }
    go(cell)*/

    //bw.close()


    /*val R = RClient() // initialise an R interpreter
    R.evalD2("as.matrix(eurodist)")
    R.x = 3
    R.matrix(3,2)=cell.codeTable // send x to R

    R.eval("mod <- glm(y~x,family=poisson())") // fit the model in R
    // pull the fitted coefficents back into scala
    val beta = DenseVector[Double](R.evalD1("mod$coefficients"))
*/

    //cells
  }

  def writeToFile(filename:String, content:String): Unit ={
    val file = new File(s"C:\\Users\\feroc\\OneDrive\\Dokumente\\Thesis\\$filename.html")
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



/*object Earth extends SimpleSwingApplication {
  val ui = new BorderPanel {
    //content
  }

  def top = new MainFrame {
    title = "title"
    contents = ui
  }

  val auth = new InputParamsDialog().auth.getOrElse(throw new IllegalStateException("You should login!!!"))
}*/

import scala.swing.BorderPanel.Position._

case class Auth(parameters:List[Any])

class InputParamsDialog extends Dialog {
  var auth: Option[Auth] = None

  title = "Parameters"
  modal = true
  var aaRSnumb =  new TextField(3)
  var aaRSlength = new TextField(3)
  var aminoAcidNumb = new CheckBox(){text = "all"}
  aminoAcidNumb.selected = true
  var codonNumb = new ComboBox(List("16 codons with length 2", "64 codons with length three", "48 strong commafree codons"))



  contents = new BorderPanel {

    layout(new FlowPanel(FlowPanel.Alignment.Center)(


      new Label("Number of aaRS genes:"),
      aaRSnumb,
      new Label("Length of an aaRS gene:"),
      aaRSlength,

      new Label("Number of amino acids:"),
      aminoAcidNumb,
      new Label("Codon set:"),
      codonNumb,

      Button("start simulation") {
        if (makeLogin()) {

          auth = Some(Auth(List(aaRSnumb.text.toInt, aaRSlength.text.toInt, aminoAcidNumb.selected, codonNumb.selection.toString)))
          close()
        } else {
          Dialog.showMessage(this, "Wrong username or password!", "Login Error", Dialog.Message.Error)
        }
      }

    ))= South


  }

  def makeLogin() = true // here comes you login logic

  centerOnScreen()
  open()
}