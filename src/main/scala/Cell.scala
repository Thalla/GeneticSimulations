import AA._

class Cell (val messengerRNA:Vector[Vector[(Int,Int)]], val geneTable:Seq[Tuple2[AA.Value, (Int, Int)]], val transferRNAs:Vector[TRNA], val aaRSs:Vector[AARS], allAARS:Vector[AARS], val genID:Int){

 /* def divide(): Cell ={
    val newAARSs:Seq[AARS] = (for{
      i <- 0 to 15          //CAUTION: 15 has to be replaced with start number of aa
    } yield (aaRSs(i).translateNew(geneTable)))
    new Cell(messengerRNA, geneTable, transferRNAs, newAARSs.toVector, genID+1)
  }*/

  def translate( ): Vector[AARS] ={
    var aaSequences: Vector[Vector[AA.Value]] = Vector()
    for(gene <- messengerRNA) {
      var sequence:Vector[AA.Value] = Vector()
        for(codon <- gene) {
          val index:Int = geneTable.indexWhere(x => x._2 == codon)
          if(index >= 0) sequence = sequence :+ geneTable(index)._1 //TODO: else case, break if aaRS is too short? Can aaRS have different lengths? Sure.
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
        case _ => (s"\nGeneration $genID: $toPrint")
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





}
