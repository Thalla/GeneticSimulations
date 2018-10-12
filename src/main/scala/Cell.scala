import Base.Base

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
class Cell (val mRNA:Vector[Vector[(Base, Base)]], val codeTable:Seq[Tuple2[AA.Value, (Base, Base)]], val transferRNAs:Vector[TRNA], val aaRSs:Vector[AARS], allAARS:Vector[AARS], val generationID:Int){
  /**
    *
    * @return
    */
  def translate( ): Vector[AARS] ={
    var aaSequences: Vector[Vector[AA.Value]] = Vector()
    for(gene <- mRNA) {
      var sequence:Vector[AA.Value] = Vector()
        for(codon <- gene) {
          val index:Int = codeTable.indexWhere(x => x._2 == codon)
          if(index >= 0) sequence = sequence :+ codeTable(index)._1 //TODO: else case, break if aaRS is too short? Can aaRS have different lengths? Sure.
      }
      aaSequences = aaSequences :+ sequence
    }
    //println(aaSequences.toString())

    //Take all aaRS that have a protein sequence fitting to one of the new aaSequences
    //val newAARSs = allAARS.filter(x => aaSequences.count(_ == x.aaSeq) >= 1)
    var newAARSs:Vector[AARS] = Vector()
    newAARSs = newAARSs ++ (allAARS find (x => aaSequences.contains(x.aaSeq)))
    println(newAARSs.toString())
    newAARSs
  }




  override def toString(): String ={
    var toPrint:String = ""
    def go(aaRSsToBePrinted: Seq[AARS]): String ={
      aaRSsToBePrinted match {
        case h :: t => toPrint +=  h.toString(); go(t)
        case _ => (s"\nGeneration $generationID: $toPrint")
      }
    }
    go(aaRSs.toList)

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
