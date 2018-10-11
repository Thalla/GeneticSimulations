import java.io.{BufferedWriter, File, FileWriter}

import org.biojava3.core.sequence.DNASequence
import AA._

import scala.util.Random

/**
  * first simple model, uses twoTuples instead of codons
  * @author Hanna Schumacher
  * @version 0.3 (previous: Simple1; less pretty code???)
  */
object Earth{

  def main(args: Array[String]): Unit ={
    //TODO: check if start parameters make sense
    //TODO: tests, exceptions
    //TODO: FIX aa for aaRS!! always use variable and not AA.Values!
    val (cell, allAARS) = init(16, 9260, 3, AA.values.toVector)  // To use each aa in the beginning start with as many aaRS as there are aa. It is not allowed to have more aaRS than aa (would be pretty useless).
    //val (cell, allAARS) = init(16, 2, Vector(Lys, Ala, Pro, Val, Thr))
    simulate(cell, allAARS)
  }


  /**
    *   creates a start cell, each cell is one generation
    *   @param aaRSnumb Same or lower as numb of possible twoTuples and lower than initAAnumb. Number of how many different aaRS proteins exist initially. Not higher than number of existing amino acids! Otherwise change code of gene table and aaRSs creation.
    *   @param aaRSlength the number of amino acids an aaRS consists of
    *   @param initAAnumb number of initially existing and used amino acids,  must be lower or same as maxTupleNumb(there are 16 twoTuple) (a tuple/codon isn't allowed to code for more than one amino acid. Any cases that do so are exspected to fail because of total chaos)
    *   @return start cell and all aaRS that can exist
    *   TODO: Add Cases (current Case0 -> each existing Tuple/Codon has one predefined aaRS): Case1 -> There can be different aaRS for one codon, leading to one codon encoding different aa leading to chaos if the aaRS exist at the same time. Case2 -> ...4*
  */
  def init (aaRSnumb:Int, possibleAARSnumb:Int, aaRSlength:Int, initAA:Vector[AA]):(Cell,Vector[AARS]) = {
    // codons
    var codons = getCodons(16)

    // mRNA
    val mRNA:Vector[Vector[(Int, Int)]] = getRandomMRNA(codons.toVector, aaRSlength, aaRSnumb).asInstanceOf[Vector[Vector[(Int, Int)]]]


    //generate code table
    var twoTuplesReduced = codons
    var initAAreduced = initAA
    val geneTable = for {
      _ <- 0 to aaRSnumb-1
    } yield {
      val rand = getUniqueRandElem(codons.seq, twoTuplesReduced.seq)
      val randAA = getUniqueRandElem(initAA, initAAreduced)
      twoTuplesReduced = twoTuplesReduced.seq.diff(Seq(rand))
      initAAreduced = initAAreduced.seq.diff(Seq(randAA)).toVector
      (randAA.asInstanceOf[AA.Value], rand.asInstanceOf[(Int,Int)])}


    // tRNA set: create a tRNA for each possible tuple
    var tRNAs:Vector[TRNA] = Vector()
     for{
      i <- 0 to (codons.length-1)
    } tRNAs = tRNAs :+ new TRNA(getUniqueRandElem(codons.seq, twoTuplesReduced).asInstanceOf[(Int, Int)], codons(i)) //TODO



    //generate tRNAs for aaRSs
    var tRNAsReduced = tRNAs
    val getRandTRNA = () => {
        //get random tRNA
        val randTRNA:TRNA = tRNAsReduced(Random.nextInt(tRNAsReduced.length))
        //delete this tRNA from tRNA list
        tRNAsReduced = tRNAsReduced.filter(_ != randTRNA)
        randTRNA
    }

    val tRNAsForCurrentAARSs:Vector[TRNA] = tRNAs.filter(x => geneTable.exists(h => (h._2 == x.anticodon)))  // choose the tRNAs that have the anticodons that fit to the current genetic Code


    // aaRS set: create a list of newly generated aaRS
    val allPossibleAaSequences = initAA.toList.combinations(aaRSlength)
    val getRandAA = () => {initAA(Random.nextInt(initAA.length))}     // Get a random aminoacid
    val getAaSeq = () => (Seq.fill(aaRSlength){getRandAA()}).toVector  // Get a random aminoacid sequence
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
    val currentAaSequences = getUniqueAaSequences(aaRSnumb)
    val additionalAaSequences = getUniqueAaSequences(possibleAARSnumb)
    var usedAA:Vector[AA] = Vector()
    def getAaAtGeneTablePos (pos:Int):AA.Value ={
      usedAA = usedAA :+ geneTable(pos)._1
      geneTable(pos)._1
    }
    //TODO: currently: goes through gene table, first aaRS has first aa and first to anticodon fitting tRNA. Desired: tRNA knows multiple anticodons!!!!
    val aaRSs:Vector[AARS] = Seq.tabulate(aaRSnumb){(counter) => {new AARS(currentAaSequences(counter), getAaAtGeneTablePos(counter), tRNAsForCurrentAARSs.filter(x => x.anticodon == geneTable(counter)._2)(0), counter)}}.toVector

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

    //TODO: currently there are as many additional aaRS produced as there are unused aa.
    // tRNAsReduced can be chosen instead of notUsedAA.length. Then the method getUniqueRandomAA has to be changed so that after each aa was chosen once a random aa is chosen
    var allAaRS:Vector[AARS] =  Seq.tabulate(possibleAARSnumb){(counter) => {new AARS(additionalAaSequences(counter), getUniqueRandomAA(), tRNAsReduced(counter%tRNAsReduced.length), aaRSnumb+counter)}}.toVector
    allAaRS = allAaRS ++ aaRSs


    println("Genetic Code:\n" + geneTable.toString())
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
   /* println("codons:\n" + codons.toString())
    println("Genetic Code:\n" + geneTable.toString())
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
    println(allAaRS.mkString(" "))*/

    var seqList:List[List[AA]] = List()
    allAaRS.foreach( x => seqList = seqList :+ x.aaSeq.toList)
    var allAaRSs = seqList.distinct


    val file = new File("C:\\Users\\feroc\\OneDrive\\Dokumente\\Thesis\\allAaRs.txt")
    val bw = new BufferedWriter(new FileWriter(file))
    bw.write(allAaRS.toString())
    bw.close()

     (new Cell(mRNA, geneTable, tRNAs, aaRSs, allAaRS, 0), allAaRS)
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
          val newGeneTable:Seq[(AA.Value, (Int, Int))] = Seq()
            (cell.geneTable.foreach(row => {
            var index:Int = newAARSs.indexWhere(aars => row._2 == aars.tRNA.anticodon)
            if(index >= 0) {
              var newRow: Seq[(AA.Value, (Int, Int))] = Seq((newAARSs(index).aa, row._2))
              newGeneTable :+ newRow
            }

          }))

          if(newGeneTable == cell.geneTable) println("Goal reached" + newGeneTable.toString())
          if((generationID+1)%50 == 0) println(generationID + "  " +newGeneTable.toString())
          var nextCell:Cell = new Cell(cell.messengerRNA, newGeneTable, cell.transferRNAs, newAARSs, allAARS, generationID+1)

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



  // codons: generate from nucleobases, codonNumb 16-> twoTuples, 64 -> codons, 48 -> scf codons
  def getCodons(codonNumb:Int):IndexedSeq[(Any)] = {
    val range:Range = codonNumb  match{
      case 16 => 0 to 0
      case 64 | 48 => 0 to 3
    }
    val filterScf = (k:Int, i:Int) => {
      if(codonNumb == 48) k!=i
      else true
    }
    var codons = (for {
      i <- 0 to 3 // four nucleobases
      j <- 0 to 3
      k <- range
      if filterScf(k, i)  // if scf option 48 then filter codons
    } yield {
      codonNumb match{
        case 16 => (i,j)
        case 64 | 48=> (i, j, k)
      }})
    codons
  }


  def getUniqueRandElem(n:Seq[Any], nReduced:Seq[Any]):Any ={
    if (nReduced.length >= 1) {
      //get random Tuple
      nReduced(Random.nextInt(nReduced.length))
    }
    else { //all elements have been used once, so they can now be used more often
      n(Random.nextInt(n.length))
    }
  }


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




  //var cells = simulate(cell)

  //println(toString(cells))






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