import akka.actor._

object EncodedAaTracker {
  def props(encodedAA:Array[Array[Int]]) : Props = Props(new EncodedAaTracker(encodedAA))
  final case class Update(aaPos:Int, hasTranslation:Boolean)
}


class EncodedAaTracker (encodedAA:Array[Array[Int]]) extends Actor {
  import EncodedAaTracker._
  var comparisons:Array[Array[Boolean]] = Array.ofDim[Boolean](3,20) //TODO
  var encodedAaChanges:Array[Array[Int]] = Array.ofDim[Int](1000000,20) //TODO
  comparisons(0) = Array.fill(20){true}
  comparisons(0)(19) = false
  var oldGen = 0
  var newGen = 1
  def receive = {
    case Update(aaPos, hasTranslation) => comparisons(newGen)(aaPos) = hasTranslation
    case genID:Int => {
      encodedAaChanges(genID-1) = Array.fill(20){0}
      for(
        i <- 0 until 19
      ){
        (comparisons(oldGen)(i), comparisons(newGen)(i)) match {
          case (false, false) => encodedAaChanges(genID-1)(i) = 0
          case (false, true) => encodedAaChanges(genID-1)(i) = 1
          case (true, false) => encodedAaChanges(genID-1)(i) = -1
          case (true, true) => encodedAaChanges(genID-1)(i) = 0
        }
      }
    }
    case Array() => encodedAaChanges
  }
}