import AA.AA
import Base.Base
import Earth.{print2DimTRNAmatrix, print3DimAARSmatrix, printMRNA}
import PrintElem.PrintElem

import scala.util.Random

/** Currently Cell doesn't aim to simulate a real cell (no cell division) but holds the state of the simulation.
  * Each cell is a new state with no connection to the old state.
  * Cell is responsible for computing the next simulation state by translating mRNA, getting a new set of aaRS and creating a new cell with it
  * TODO connect states
  * @param mRNA to be translated to get the aaRS and therefore the codeTable for the next cell generation
  * @param codeTable either check the code table or the current aaRSs themselves to find out about the translation from codon to amino acid
  * @param transferRNAs not relevant yet
  * @param aaRSs current set of aaRS, used for translating the mRNA to get the aaRS set for the next cell generation
  * @param allAARS this list holds all valid aaRS that could exist and their randomly predefined behaviour
  * @param generationID
  */
class Cell (val mRNA:List[List[Int]], var allTRNAs: Array[Array[TRNA]], var aaRSs:Vector[Vector[AA.Value]], val allAARS:Array[Array[Array[AARS]]],  val initAA:Vector[AA], val generationID:Int){
  /**
    *
    * @return
    */
  def translate( ): Cell ={
    var newAARSs:Vector[Vector[AA.Value]] = Vector()
    var aaSequences: Vector[Vector[AA.Value]] = Vector()
    for(gene <- mRNA) {
      var sequence:Vector[AA.Value] = Vector()
        for(codon <- gene) {
           //to codon fitting tRNAs that are used by the current aaRS set
          val fittingTRNAsInUse:Array[TRNA] = allTRNAs(codon).filter(tRNA => if(tRNA == null){false}else{tRNA.aaRSsequences.length > 0})

          val tRNA = fittingTRNAsInUse(Random.nextInt(fittingTRNAsInUse.length)) //random TRNA from the fitting list
          val aaSeq:Vector[AA] = tRNA.aaRSsequences(Random.nextInt(tRNA.aaRSsequences.length))
          val aaRS:AARS = allAARS(aaSeq(0).id)(aaSeq(1).id)(aaSeq(2).id)
          sequence = sequence :+ aaRS.aa(Random.nextInt(aaRS.aa.length))              //TODO: else case, break if aaRS is too short? Can aaRS have different lengths? Sure. But not important? There are no more than 20^3 functional aaRS
      }
      var newAARS:AARS = allAARS(sequence(0).id)(sequence(1).id)(sequence(2).id)
      if(newAARS != null){
        if(aaRSs.contains(newAARS)) {  //lifeticks are not 0, remove aaRS from old list to not have it twice afterwards and to not reduce the lifeticks
          aaRSs = aaRSs.filter(_ != newAARS.aaSeq)
        }
        else{
          //activate tRNAs
          newAARS.tRNAs.foreach(tRNA => allTRNAs(tRNA._1)(tRNA._2).addAARSsequence(newAARS.aaSeq))
        }
        newAARS.resetLifeticks()
        newAARSs = newAARSs :+ newAARS.aaSeq
      }
      else{
        var aminoAcids:Vector[AA] = Vector()
        var tRNAs:Vector[(Int,Int)] = Vector()

        //AA
        val numberOfAA = Random.nextInt(2)
        for(
          i <- 0 to numberOfAA
        ){
          aminoAcids = aminoAcids :+ initAA(Random.nextInt(initAA.length))
        }

        //tRNA
        val numberOfTRNAs = Random.nextInt(2)
          for(
            i <- 0 to numberOfTRNAs
          ){
          val anticodonPos = Random.nextInt(allTRNAs.length)
          val stemPos = Random.nextInt(allTRNAs.length)
          if(allTRNAs(anticodonPos)(stemPos) == null){
            allTRNAs(anticodonPos)(stemPos) = new TRNA(anticodonPos, stemPos, Vector(sequence))
          }
            tRNAs = tRNAs :+ (anticodonPos, stemPos)
        }

        newAARS = new AARS(sequence,tRNAs,aminoAcids)

        allAARS(sequence(0).id)(sequence(1).id)(sequence(2).id) = newAARS
        newAARSs = newAARSs :+ newAARS.aaSeq
      }
    }
    aaRSs.foreach(seq => {
      allAARS(seq(0).id)(seq(1).id)(seq(2).id).reduceLifeTicks()
      if(allAARS(seq(0).id)(seq(1).id)(seq(2).id).lifeticks <= 0){ //if aaRS is dead
        allAARS(seq(0).id)(seq(1).id)(seq(2).id).tRNAs.foreach(tRNA => allTRNAs(tRNA._1)(tRNA._2).aaRSsequences.filter(_!=seq))  //delete tRNA connection to this aaRS
      }
    })

    newAARSs = newAARSs ++ aaRSs

    //print initial stuff
    def printElements(toPrint:List[PrintElem]):Unit = {
      toPrint.foreach(elem => {
        elem match {
          case PrintElem.allTRNA => print2DimTRNAmatrix(allTRNAs)
          case PrintElem.aaRSs => aaRSs.toString()
          case PrintElem.allAARS => print3DimAARSmatrix(allAARS)
        }
      })
    }
    printElements(List(PrintElem.allTRNA, PrintElem.aaRSs, PrintElem.allAARS))
    new Cell(mRNA, allTRNAs, newAARSs, allAARS, initAA, generationID+1)
  }

  /*def reduceLifeTicks():Vector[AARS]={
    var updatedAARS:Vector[AARS] = for(
      aaRS<- aaRSs
    )yield{
      //work with copy
      aaRS.reduceLifeTicks()
      if(aaRS.lifeticks <= 0){
        allTRNAs = aaRS.apoptosis(allTRNAs)
      }
      //make copy real, replace old version
      allAARS(aaRS.aaSeq(0).id)(aaRS.aaSeq(1).id)(aaRS.aaSeq(2).id) = aaRS
      aaRS
    }
    updatedAARS.filter(aaRS => aaRS.lifeticks > 0)
  }*/






  override def toString(): String ={
    var toPrint:String = ""
    def go(aaRSsToBePrinted: Seq[AARS]): String ={
      aaRSsToBePrinted match {
        case h :: t => toPrint +=  h.toString(); go(t)
        case _ => (s"\nGeneration $generationID: $toPrint")
      }
    }
    go(List())

  }












  // same as above (create aaRSs)
  /*for (i <- 1 to aaRSnumb){
    var aaSeq:Vector[AA] = Vector()

    for(j <- 1 to aaRSlength) {
      aaSeq = aaSeq :+ AA(scala.util.Random.nextInt(AA.maxId))
    }
    aaRSs = aaRSs :+ new AARS(aaSeq, i)
  }*/


  //println("AARSs")
  //println((aaRSs(1)).toString())



  /* def divide(): Cell ={
     val newAARSs:Seq[AARS] = (for{
       i <- 0 to 15          //CAUTION: 15 has to be replaced with start number of aa
     } yield (aaRSs(i).translateNew(geneTable)))
     new Cell(messengerRNA, geneTable, transferRNAs, newAARSs.toVector, genID+1)
   }*/




}
