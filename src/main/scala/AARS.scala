import AA._

class AARS(val aaSeq: Vector[AA], val tRNAs:Vector[(Int,Int)], val aa:Vector[AA]) {

  var lifeticks = 2  //be careful to also change resetLifeticks if initial value is changed

  def reduceLifeTicks():Unit ={
    lifeticks -= 1
  }

  def resetLifeticks():Unit ={
    lifeticks = 2
  }

  // when an aaRS dies the usage of its tRNAs stops, they are unloaded and can't recognize any codons
  def apoptosis(allTRNA:Array[Array[TRNA]]):Array[Array[TRNA]]={
    tRNAs.foreach(tRNA => {
      allTRNA(tRNA._1)(tRNA._2).removeAARSsequence(aaSeq)
    })
    allTRNA
  }




  override def toString(): String ={
    "\n" + "aaRS: " + aaSeq.mkString(", ") + "   " + aaSeq.toString() +  "TRNAs: " + tRNAs.mkString(", ")
  }



}

