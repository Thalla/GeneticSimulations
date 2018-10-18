
import java.io.{BufferedWriter, File, FileWriter}

import PrintElem.PrintElem

import scala.swing._
import scala.swing.event.ButtonClicked
//import org.biojava3.core.sequence.DNASequence
import AA._
import Base.Base
import scala.util.Random


/** first simple model
  * uses twoTuples instead of codons, doesn't make use of tRNA stem
  * made for 16 initial aaRS (for each twoTuple one) with random amino acid sequences of length 3, 20 amino acids to choose from, nearly every possible amino acid sequence is a valid aaRS
  * @problem The next generation can have an aaRS recognizing anticodon XY and binding amino acid Abc as well as an aaRS also recognizing XY but binding it to Hij.
  *          So there are two translations for one codon and the current code table holds just one
  * @solutions V1: The case is not allowed because it would lead to chaos.
  *            V2: The case is allowed because it is allowed in nature (TODO: some more research) .
  *                => The code table must hold multiple translations.
  *            Either way the code table must be extended to :TODO: include the last generations
  *           (Delete code table and just use AARS? -> make AARS easier to use)
  *           TODO: Define clear relationship between AARS and tRNA (n:m?)
  * @codeGoal flexible NN construct?
  * @author Hanna Schumacher
  * @version 0.5 : more documentation (descriptions, undefined cases, param impact), added Base Enumeration (previous: GeneticSimulations)
  */
object Earth{

  //TODO: check if start parameters are valid
  //TODO: tests, exceptions, logging, class diagram, cases, directly connect code/model decisions to research
  //TODO: FIX aa for aaRS!! always use variable and not AA.Values! (Check in Cell too)

  /** initiates and starts simulation
    * @param args
    */
  def main(args: Array[String]): Unit ={

    //val auth = new InputParamsDialog().auth.getOrElse(throw new IllegalStateException("You should login!!!"))
    //println(auth.toString())
    simulate(init())           //
  }


  /** creates a start cell, each cell is one generation
    * @param aaRSnumb Number of how many different aaRS proteins exist initially.
    *                 Must be same or lower as numb of possible twoTuples because aaRSnumb defines the length of the code table. Therefore aaRSnumb isn't allowed to be bigger than max numb of codons/tuples.
    *                 It must be lower than initAA.length/nunmber of existing amino acids. Otherwise the same AA is coded by various aaRS.
    *                 aaRSnumb defines the length of the code table.  TODO: allow codons to be translated by multiple aaRS (V1: to multiple AA (-> chaos?), V2: to always the same AA) -> changes in gene table and aaRSs creation
    * @param aaRSlength The number of amino acids an aaRS consists of. TODO: aaRS can have varying lengths
    * @param initAA initially existing and used amino acids, must be lower or same as maxTupleNumb(there are 16 twoTuple) TODO: allow having new/non proteinogenic amino acids
    * @param codonNumb Either 16 (set of twoTuples is created) or 64 (set of real codons is created) or 48 (set of strong commafree codons is created)
    *  @return start cell and all aaRS that can exist
  */
  def init (aaRSnumb:Int = 20, aaRSlength:Int = 3, initAA:Vector[AA] = AA.values.toVector, codonNumb:Int = 16):Cell= {
    // create codons
    val codons = getCodons(codonNumb)

    // start connection is: codon with id x has tRNA with id x has aaRS with id x (three tRNA have a second aaRS) has amino acid with id x

    // create random mRNA
    val mRNA:List[List[Int]] = getRandomMRNA(codons.toList, aaRSlength, aaRSnumb)


    //create tRNAs: for each codon one (no codon anticodon translation, both are same)
    var allTRNA =  Array.ofDim[TRNA](codons.length, codons.length)  //TODO needed???
    var tRNAs:Vector[(Int, Int)] = Vector.tabulate(codons.length)(i => {
      val stemCodonPos = Random.nextInt(codons.length)
      val tRNA =  new TRNA(i, stemCodonPos, Vector()) //leere Hülle
      allTRNA (i)(stemCodonPos) = tRNA
      (i, stemCodonPos)
    })
    for(
      i <- 0 until aaRSnumb - codons.length
    ){
      val anticodonPos = Random.nextInt(codons.length)
      val stemCodonPos = Random.nextInt(codons.length)
      var tRNA:TRNA =  new TRNA(anticodonPos, stemCodonPos, Vector())
      tRNAs = tRNAs :+ (anticodonPos, stemCodonPos) //TODO Fehler wenn gleiche TRNA wie schon erzeugt
      allTRNA (anticodonPos)(stemCodonPos) = tRNA
    }

    //create aaRSs
    var allAARS =  Array.ofDim[AARS](20,20, 20) //TODO ? more dimensions ?
    var aaRSs:Vector[Vector[AA.Value]] = Vector()
    //create amino acid sequences
    val aaRSsequences:IndexedSeq[IndexedSeq[AA]] = for(
      i <- 0 until aaRSnumb
    )yield{
      val seq:IndexedSeq[AA] = for(
        j <- 0 until aaRSlength
      ) yield AA(Random.nextInt(initAA.length))
      seq
    }
    // create an aaRS from each amino acid sequence in aaRSsequences and store it in an aaRS Vector (for current cell) and Matrix (collects all aaRS)
    var counter:Int = 0
    aaRSsequences.foreach(aaSeq => {
      var x = tRNAs(17)
      var y = allTRNA(tRNAs(17)._1)(tRNAs(17)._2)
      allTRNA(tRNAs(counter)._1)(tRNAs(counter)._2).addAARSsequence(aaSeq.toVector)
      val aaRS = new AARS(aaSeq.toVector, tRNAs, Vector(initAA(counter))) //TODO aaRSnumb > initAA.length
      aaRSs = aaRSs :+ aaRS.aaSeq
      allAARS(aaSeq(0).id) (aaSeq(1).id) (aaSeq(2).id) = aaRS
      counter += 1
    })


    //print initial stuff
    def printElements(toPrint:List[PrintElem]):Unit = {
      toPrint.foreach(elem => {
        elem match {
          case PrintElem.codons => println("Codons" + codons.toString())
          case PrintElem.mRNA => printMRNA(mRNA, codons)
          case PrintElem.tRNAs => println("tRNAs:\n" + tRNAs.toString()); println()
          case PrintElem.allTRNA => print2DimTRNAmatrix(allTRNA)
          case PrintElem.aaRSs => aaRSs.toString()
          case PrintElem.allAARS => print3DimAARSmatrix(allAARS)
        }
      })
    }
    printElements(List(PrintElem.codons, PrintElem.mRNA, PrintElem.tRNAs, PrintElem.aaRSs))

    def getHTMLstring(toHTML:List[PrintElem]):String = {
      var content = ""
      toHTML.foreach(elem => {
        elem match {
          case PrintElem.codons => {
            content += "<h2>Codons</h2>\n"
            content += codons.mkString("<p>",", ","</p>") + "</br>\n"
          }
          case PrintElem.mRNA => content += "<h2>mRNA</h2>\n"
            content += mRNAtoHTML(mRNA, codons)
          case PrintElem.tRNAs => println("tRNAs:\n" + tRNAs.toString()); println()
          case PrintElem.allTRNA => print2DimTRNAmatrix(allTRNA)
          case PrintElem.aaRSs => aaRSs.toString()
          case PrintElem.allAARS => print3DimAARSmatrix(allAARS)
        }
      })
      content
    }
    val content = getHTMLstring(List(PrintElem.codons, PrintElem.mRNA, PrintElem.tRNAs, PrintElem.aaRSs))




    val file = new File("C:\\Users\\feroc\\OneDrive\\Dokumente\\Thesis\\webpage.txt")
    val bw = new BufferedWriter(new FileWriter(file))
    bw.write(content.toString())
    bw.close()

    new Cell(mRNA, allTRNA, aaRSs, allAARS, initAA, 0)

  }




  // simulation; holds global parameters; cares for cell reproduction
  def simulate(cell:Cell): Vector[Cell] = {
    //Each cell represents one generation and Earth holds all.
    var cells: Vector[Cell] = Vector()
    cells = cells :+ cell

    def go(cell:Cell):Unit = {
      var newCell:Cell = cell.translate()
      cells = cells :+ newCell
      go(newCell)
    }
    go(cell)

    cells
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

  def print2DimTRNAmatrix (matrix:Array[Array[TRNA]]):Unit = {
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



//aaRSs = aaRSs :+ aaRSs(0)
//              //dup.groupBy(identity).collect { case (x,ys) if ys.lengthCompare(1) > 0 => x }
//    aaRSs = (aaRSs.groupBy(identity).collect { case (x, Vector(_,_*)) => x }).toVector




/*
// twoTuples: generate from nucleobases
val twoTuples: IndexedSeq[(Int, Int)] = (for {
  i <- 0 to 3 // four nucleobases
  j <- 0 to 3
} yield (i, j))
println("twoTuples:\n" + twoTuples.toString())

// codons: generate from nucleobases
val codons: IndexedSeq[(Int, Int, Int)] = (for {
  i <- 0 to 3 // four nucleobases
  j <- 0 to 3
  k <- 0 to 3 if(k != i)
} yield (i, j, k))
println("codons:\n" + codons.toString())*/

/*    //generate tRNAs for aaRSs
    var tRNAsReduced = tRNAs
    val getRandTRNA = () => {
        //get random tRNA
        val randTRNA:TRNA = tRNAsReduced(Random.nextInt(tRNAsReduced.length))
        //delete this tRNA from tRNA list
        tRNAsReduced = tRNAsReduced.filter(_ != randTRNA)
        randTRNA
    }*/



// funktioniert:
/*val getRandAA = () => {initAA(Random.nextInt(initAA.length))}     // Get a random aminoacid
val getAaSeq = () => (Seq.fill(aaRSlength){getRandAA()}).toVector // Get a random aminoacid sequence
val getUniqueAaSequences = (numb:Int) => {
      var aaSequences: Vector[Vector[AA.Value]] = Vector()
      @annotation.tailrec
      def go(counter: Int) :Vector[Vector[AA.Value]] = {
        var aaSeq = getAaSeq()
        //prevents having the same sequence twice
        while(aaSequences.contains(aaSeq)){
          aaSeq = getAaSeq()
        }
        // adds sequence to list
        aaSequences = aaSequences :+ aaSeq
        // checks if desired number of sequences was generated
        if(counter >= numb){
          aaSequences
        }
        else{
          go(counter+1)
        }
      }
      go(0)
    }
    val currentAaSequences = getUniqueAaSequences(aaRSnumb)*/


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