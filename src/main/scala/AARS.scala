import AA._

class AARS(val aaSeq: Vector[AA],  val aa:AA, val tRNA:TRNA, id: Int) {

  // Takes: mRNA for aaRS production
  // Does: reads mRNA and translates codon by codon into amino acids
  // Uses: genetic Code as basic for translation
  // Result:
  def translateNew(geneticCode:Seq[Tuple2[AA.Value, (Int, Int)]]): AARS ={

    //var mRNAtoTranslate = mRNA
    var newAaSeq:Vector[AA]= Vector()

    def go (newAaSeq:Vector[Any]): Vector[AA] ={
      //var tuple:(Int, Int) = mRNAtoTranslate.head
      //mRNAtoTranslate = mRNAtoTranslate.tail
      geneticCode match{
        case Seq(a,(tuple)) => go(newAaSeq :+ a)
      }
    }
    go(newAaSeq)

    def calculateNewAA(tRNA:TRNA, newAaSeq:Vector[AA]): AA ={
      //TODO
      Lys
    }

    new AARS(newAaSeq, aa, tRNA, id)

  }


  override def toString(): String ={
    "\n" + id + ". aaRS: " + aaSeq.mkString(", ") + ",    bindingAA: " + aa + ",      bindingTRNA:  " + tRNA.stem + tRNA.anticodon
  }



}

