import AA._

import scala.util.Random
import java.io.{BufferedWriter, File, FileWriter}

import com.github.tototoshi.csv._

import scala.collection.mutable.{ArrayBuffer, ListBuffer}  //TODO?: replace by own methods (data-> List-> needed type => data-> needed type)


/** Fourth draft
  *
  * @nextSteps: similarity, probabilities, NN?, DM?, evo. Alg.?, analyse, parameters, Empirie, Fitness, Metrik
  * @author Hanna Schumacher
  * @version 4.1   -> 500.000 in 315s ; in 2?? s; in 272s without some stuff, 190s without all
  * TODO choose by translation fitness, not always the best one, use pathvariable
  */
object Earth{
  /*
  def writeToFile(filename:String, content:String): Unit ={
    val file = new File(s"C:\\Users\\feroc\\OneDrive\\Dokumente\\HS\\Semester\\4\\Thesis\\$filename.html")
    val bw = new BufferedWriter(new FileWriter(file))
    bw.write(content)
    bw.close()
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



  def print2DimTRNAmatrix (matrix:Array[Array[Any]]):Unit = {
    for (
      i <- 0 until matrix.length;
      j <- 0 until matrix(0).length
    ) {
      val elem = matrix(i)(j)
      if (elem == null) print(" ")//i + ", " + j + " ")
      else println(matrix(i)(j).toString())
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
*/
}