import AA._

class TRNA (val anticodon: Int, val stem:Int,var aaRSsequences:Vector[Vector[AA]]) {


  def addAARSsequence(seq:Vector[AA]):Unit={
    aaRSsequences = aaRSsequences :+ seq
  }


  def removeAARSsequence(seq:Vector[AA]):Unit={
    aaRSsequences = aaRSsequences.filter(_ != seq)
  }


  override def toString(): String = {
    "TRNA: " + stem.toString() +" "+ anticodon.toString()
  }
}
