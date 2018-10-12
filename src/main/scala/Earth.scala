import java.io.{BufferedWriter, File, FileWriter}
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
  * @version 0.4 : more documentation (descriptions, undefined cases, param impact), added Base Enumeration (previous: Simple2)
  */
object Earth{

  //TODO: check if start parameters are valid
  //TODO: tests, exceptions, logging, class diagram, cases, directly connect code/model decisions to research
  //TODO: FIX aa for aaRS!! always use variable and not AA.Values! (Check in Cell too)

  /** initiates and starts simulation
    * @param args
    */
  def main(args: Array[String]): Unit ={
    val (cell, allAARS) = init(16, 9260, 3, AA.values.toVector)  // To use each aa in the beginning start with as many aaRS as there are aa. It is not allowed to have more aaRS than aa (would be pretty useless).
    //val (cell, allAARS) = init(16, 2, Vector(Lys, Ala, Pro, Val, Thr))
    simulate(cell, allAARS)
  }


  /** creates a start cell, each cell is one generation
    * @param aaRSnumb Number of how many different aaRS proteins exist initially.
    *                 Must be same or lower as numb of possible twoTuples because aaRSnumb defines the length of the code table. Therefore aaRSnumb isn't allowed to be bigger than max numb of codons/tuples.
    *                 It must be lower than initAA.length/nunmber of existing amino acids. Otherwise the same AA is coded by various aaRS.
    *                 aaRSnumb defines the length of the code table.  TODO: allow codons to be translated by multiple aaRS (V1: to multiple AA (-> chaos?), V2: to always the same AA) -> changes in gene table and aaRSs creation
    * @param aaRSlength The number of amino acids an aaRS consists of. TODO: aaRS can have varying lengths
    * @param initAA initially existing and used amino acids, must be lower or same as maxTupleNumb(there are 16 twoTuple) TODO: allow having new/non proteinogenic amino acids
    * @return start cell and all aaRS that can exist
  */
  def init (aaRSnumb:Int, possibleAARSnumb:Int, aaRSlength:Int, initAA:Vector[AA]):(Cell,Vector[AARS]) = {
    // create codons
    val codons = getCodons(16)

    // create random mRNA
    val mRNA:Vector[Vector[(Base, Base)]] = getRandomMRNA(codons.toVector, aaRSlength, aaRSnumb).asInstanceOf[Vector[Vector[(Base, Base)]]]

    // create code table
    // as long as aaRSnumb isn't bigger than codons.length and initAA.length, no codon or amino acid is used twice
    var codonsReduced = codons
    var initAaReduced = initAA
    val codeTable = for {
      _ <- 0 until aaRSnumb
    } yield {
      // get random elements
      val randCodon = getUniqueRandElem(codons.seq, codonsReduced.seq)
      val randAA = getUniqueRandElem(initAA, initAaReduced)
      // assure that the chosen elements won't be used in the next iteration
      codonsReduced = codonsReduced.seq.diff(Seq(randCodon))
      initAaReduced = initAaReduced.seq.diff(Seq(randAA)).toVector
      // return tuple (amino acid, (base, base))
      (randAA.asInstanceOf[AA.Value], randCodon.asInstanceOf[(Base, Base)])
    }


    // create tRNAs, one for each possible codon/twoTuple; stem is chosen randomly but no stem twice
    var tRNAs:Vector[TRNA] = Vector()
    val usedCodonPositions:Seq[(Int, Int)] = for{
      i <- 0 until codons.length
    } yield {
      //chooses randomly (from the codons used in the code table) which codons are taken for tRNA stem
      val stemPos = Random.nextInt(codons.length-1)
      tRNAs = tRNAs :+ new TRNA(codons(stemPos), codons(i))
      (stemPos, i) //collects the indices of the already used codon tuples to not use them again
      }
    // create a list with all possible tRNAs  TODO: Is it needed that tRNAs with same stem and anticodon are implemented as the same object? If not implement this easier.
    var allTRNAs: Vector[TRNA] = Vector()
    for{
      i <- 0 to (codons.length-1)
      j <- 0 to (codons.length-1)
    } {
      if(!usedCodonPositions.contains((j,i))){
        allTRNAs = allTRNAs :+ new TRNA(codons(j), codons(i))
      }
    }
     allTRNAs = allTRNAs ++ tRNAs

    var tRNAaaRS:Map[TRNA, Vector[AARS]] = Map()
    for{
      t <- 0 until tRNAs.length
    }{
      // 16 different stems, 1 codon, 16 different amino acids   ; less amino acids than codons :  cta ctaa ctaaa cta ctaa, cctaaa ccta cctaa cctaaa ccta   ; less codons than amino acids: some codons can't
      tRNAaaRS += (tRNAs(t) -> Vector(initAA(t % (initAA.length))))
    }

    // choose the tRNAs that have the anticodons that fit to the current genetic Code (for cases in which the number of initial aaRS is lower than number of possible codons TODO: don't make a tRNA (for tRNAs) for each codon  but for each codon in code table so that this step can be skipped
    val tRNAsForCurrentAARSs:Vector[TRNA] = tRNAs.filter(x => codeTable.exists(h => (h._2 == x.anticodon)))


    // create all possible aaRS and take randomly $aaRSnumb aaRSs           TODO solution for aaRS with sequence length > 3
    val allPossibleAaSequences = initAA.toList.combinations(aaRSlength)
    var currentAaSequences: Vector[Vector[AA]] = Vector()
    var additionalAaSequences: Vector[Vector[AA]] = allPossibleAaSequences.toVector.asInstanceOf[Vector[Vector[AA]]]
    val allAaSequences = additionalAaSequences
    for{
      _ <- 0 until aaRSnumb
    }{
      val randAaSeq:Vector[AA] = getUniqueRandElem(allAaSequences, additionalAaSequences).asInstanceOf[Vector[AA]]
      currentAaSequences = currentAaSequences :+ randAaSeq
      additionalAaSequences = additionalAaSequences.filter(_ != randAaSeq)
    }


    //val additionalAaSequences = getUniqueAaSequences(possibleAARSnumb)
    var usedAA:Vector[AA] = Vector()
    def getAaAtGeneTablePos (pos:Int):AA.Value ={
      usedAA = usedAA :+ codeTable(pos)._1
      codeTable(pos)._1
    }
    //TODO: currently: goes through gene table, first aaRS has first aa and first to anticodon fitting tRNA. Desired: tRNA knows multiple anticodons!
    val aaRSs:Vector[AARS] = Seq.tabulate(aaRSnumb){(counter) => {new AARS(currentAaSequences(counter), getAaAtGeneTablePos(counter), tRNAsForCurrentAARSs.filter(x => x.anticodon == codeTable(counter)._2)(0), counter)}}.toVector

    var aa = AA.Leu
    var notUsedAA:Vector[AA] = AA.values.toVector.filterNot(x => usedAA.contains(x))
    val getUniqueRandomAA = () => {

      if(notUsedAA.length > 0) {
        aa = notUsedAA(Random.nextInt(notUsedAA.length))
        notUsedAA = notUsedAA.filterNot(x => x == aa)
      }
      else{
        aa = AA(Random.nextInt(AA.maxId))
      }
      aa
    }

    //TODO: how many different tRNAs? How many stems for one codon? To how many aaRS does one tRNA fit? => see cases
    var tRNAsReduced = tRNAs
    tRNAsReduced = tRNAsReduced.filter(!tRNAsForCurrentAARSs.contains(_))
    var allAaRS:Vector[AARS] =  Seq.tabulate(possibleAARSnumb){(counter) => {new AARS(additionalAaSequences(counter), getUniqueRandomAA(), tRNAsReduced(counter%tRNAsReduced.length), aaRSnumb+counter)}}.toVector
    allAaRS = allAaRS ++ aaRSs


    println("Genetic Code:\n" + codeTable.toString())
    println("mRNA:\n" + mRNA.toString())
    println("tRNAs")
    for{
      i <- 0 to (codons.length-1)
    } print(tRNAs(i).stem.toString + tRNAs(i).anticodon.toString + "; ")

    println()
    println("tRNAsForCurrentAARSs:")
    for{
      i <- 0 to (tRNAsForCurrentAARSs.length-1)
    } print(tRNAsForCurrentAARSs(i).stem.toString + tRNAsForCurrentAARSs(i).anticodon.toString + "; ")
    println()
    println(aaRSs.mkString(" "))
   /* println("codons/twoTuples:\n" + codons.toString())
    println("random code table:\n" + geneTable.toString())
    println("random mRNA:\n" + mRNA.toString())
    println("random tRNAs")
    for{
      i <- 0 to (codons.length-1)
    } print(tRNAs(i).stem.toString + tRNAs(i).anticodon.toString + "; ")

    println()
    println("tRNAsForCurrentAARSs:")
    for{
      i <- 0 to (tRNAsForCurrentAARSs.length-1)
    } print(tRNAsForCurrentAARSs(i).stem.toString + tRNAsForCurrentAARSs(i).anticodon.toString + "; ")
    println()
    println(aaRSs.mkString(" "))
    println(allAaRS.mkString(" "))*/

    var seqList:List[List[AA]] = List()
    allAaRS.foreach( x => seqList = seqList :+ x.aaSeq.toList)
    var allAaRSs = seqList.distinct


    /*
    val file = new File("C:\\Users\\feroc\\OneDrive\\Dokumente\\Thesis\\allAaRs.txt")
    val bw = new BufferedWriter(new FileWriter(file))
    bw.write(allAaRS.toString())
    bw.close()
    */

     (new Cell(mRNA, codeTable, tRNAs, aaRSs, allAaRS, 0), allAaRS)
  }




  // simulation; holds global parameters; cares for cell reproduction
  def simulate(cell:Cell, allAARS:Vector[AARS]): Vector[Cell] = {
    //Each cell represents one generation and Earth holds all.
    var generationID: Int = 0
    var cells: Vector[Cell] = Vector()
    cells = cells :+ cell
    var fail:Int = 0


      def go(cell:Cell, aaRSs:Vector[AARS]):Unit = {
        val newAARSs = cell.translate()
        if(newAARSs.length != 0) {
          //var newAllAARS = allAARS ++ newAARSs //TODO: lazy AARS
          val newGeneTable:Seq[(AA.Value, (Base, Base))] = Seq()
            (cell.codeTable.foreach(row => {
            var index:Int = newAARSs.indexWhere(aars => row._2 == aars.tRNA.anticodon)
            if(index >= 0) {
              var newRow: Seq[(AA.Value, (Base, Base))] = Seq((newAARSs(index).aa, row._2))
              newGeneTable :+ newRow
            }

          }))

          if(newGeneTable == cell.codeTable) println("Goal reached" + newGeneTable.toString())
          if((generationID+1)%50 == 0) println(generationID + "  " +newGeneTable.toString())
          var nextCell:Cell = new Cell(cell.mRNA, newGeneTable, cell.transferRNAs, newAARSs, allAARS, generationID+1)

          cells = cells :+ nextCell
          //println("Living Generation " + generationID)

          generationID = generationID + 1
          go(nextCell, allAARS)
        }
        else{
          val (nextCell, allAARS) = init(20, 5000, 3, AA.values.toVector)//Vector(Lys, Ala, Pro, Val, Thr))//AA.values.toVector)
          generationID = 0
          fail = fail +1
          if((fail+1)%50==0) println("Fail: " + fail)
          go(nextCell, allAARS)
        }
      }
      go(cells(generationID), allAARS)



    cells
  }


  /** generate from nucleobases
    * @param codonNumb 16 -> twoTuples, 64 -> codons, 48 -> strong comma free codons
    * @return either twoTuples or codons or strong comma free codons
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
  def getRandomMRNA(codons:Vector[Any], geneLength:Int, geneNumb:Int):Vector[Vector[Any]] = {
    var mRNA:Vector[Vector[Any]] = Vector()

    val getRandomGene = () => {
      var gene:Vector[Any] = Vector()
      for {
        _ <- 0 to (geneLength-1)
      } gene = gene :+ codons(Random.nextInt(codons.length))
      gene
    }
    // get a gene that isn't already existing
    val getUniqueRandomGene = () => {
      var gene:Vector[Any] = getRandomGene()
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