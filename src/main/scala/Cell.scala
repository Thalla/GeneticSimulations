import java.io.{BufferedWriter, File, FileWriter}

import AA.AA
import Base.Base
import Earth.{mRNAtoHTML, print2DimTRNAmatrix, print3DimAARSmatrix, printMRNA, getCodons}
import PrintElem.PrintElem

import scala.util.Random

/** Cell holds the state of the simulation.
  * Each cell is a new state with no connection to the old state.
  * Cell is responsible for computing the next simulation state by translating mRNA, getting a new set of aaRS and creating a new cell with it
  * TODO connect states, reduce constructor parameters
  * @param mRNA to be translated to get the aaRS and therefore the codeTable for the next cell generation
  * @param codeTable either check the code table or the current aaRSs themselves to find out about the translation from codon to amino acid
  * @param livingAARSs current set of aaRS, used for translating the mRNA to get the aaRS set for the next cell generation
  * @param allAARS this list holds all valid aaRS that could exist and their randomly predefined behaviour
  * @param generationID nextCell has generationID+1
  */
class Cell (val mRNA:List[List[Int]], var codeTable:Array[Array[List[AARS]]], var livingAARSs:Vector[AARS], val allAARS:Array[Array[Array[AARS]]], val initAA:Vector[AA], val codonNumb:Int, val generationID:Int) {
  /**
    *
    * @return
    */
  def translate(): Cell = {

    livingAARSs.foreach(aaRS => {
      aaRS.reduceLifeTicks()
      if(aaRS.lifeticks <= 0) livingAARSs = livingAARSs.filter(_ != aaRS)
    })

    for (
      gene <- mRNA
    ) {
      var sequence: Vector[AA.Value] = Vector()
      var isStopped:Boolean = false
      for (
        codon <- gene if !isStopped
      ) {
        val usableAARS = codeTable(codon).filter(_ != null) //find all aaRS that can translate this codon
        if (!usableAARS.isEmpty) {
          val aarsList: List[AARS] = usableAARS(Random.nextInt(usableAARS.length))
          val aaRS: AARS = aarsList(Random.nextInt(aarsList.length))
          val translationKeys = aaRS.translations.keySet
          val toCodonFittingTranslationKeys = translationKeys.filter(_._2 == codon)
          val aa: AA = toCodonFittingTranslationKeys.toList(Random.nextInt(toCodonFittingTranslationKeys.size))._1
          sequence = sequence :+ aa
        }
        else {
          isStopped = true
        }
      }
      if(!isStopped){
        val newAARS = getAARS(sequence)
        if (!livingAARSs.contains(newAARS)) {
          livingAARSs = livingAARSs :+ newAARS
      }
        sequence = Vector()
      }
    }

    val newCodeTable:Array[Array[List[AARS]]] = Array.ofDim[List[AARS]](codonNumb, initAA.length)
    livingAARSs.foreach(aaRS => {
      aaRS.translations.foreach(translation => {
        if(newCodeTable(translation._1._2)(translation._1._1.id) == null) newCodeTable(translation._1._2)(translation._1._1.id)= List(aaRS)
        else newCodeTable(translation._1._2)(translation._1._1.id) = newCodeTable(translation._1._2)(translation._1._1.id) :+ aaRS
        allAARS(aaRS.aaSeq(0).id)(aaRS.aaSeq(1).id)(aaRS.aaSeq(2).id) = aaRS
      })
    })


    new Cell(mRNA, newCodeTable, livingAARSs, allAARS, initAA, codonNumb, generationID + 1)
  }


  /**
    *
    * @param sequence
    * @return
    */
  def getAARS(sequence:Vector[AA.Value]):AARS={
    var newAARS: AARS = allAARS(sequence(0).id)(sequence(1).id)(sequence(2).id)
    if (newAARS != null) {
      newAARS.resetLifeticks()
    }
    else {
      var translations:Map[(AA, Int),List[(Double, Int)]] = Map()
      val numbOfAnticodonsForOneAARS = List(1,2,3,4,5,6) //TODO remove copy paste
      val numbAnticodons= numbOfAnticodonsForOneAARS(Random.nextInt(numbOfAnticodonsForOneAARS.length))
      for(
        _ <- 0 until numbAnticodons
      ){
        translations += (initAA(Random.nextInt(initAA.length)), Random.nextInt(codonNumb))-> List((1.0,Random.nextInt(codonNumb))) //TODO choose similar codons (first random, rest +x)
      }
      newAARS = new AARS(sequence, translations)
      allAARS(sequence(0).id)(sequence(1).id)(sequence(2).id) = newAARS

    }
    newAARS
  }


  /**
    *
    * @return
    */
  override def toString(): String = {
    var toPrint: String = ""

    def go(aaRSsToBePrinted: Seq[AARS]): String = {
      aaRSsToBePrinted.toList match {
        case h :: t => toPrint += h.toString(); go(t)
        case _ => toPrint += (s"\nGeneration $generationID: $toPrint"); toPrint
      }
    }
    go(livingAARSs)

  }


  /**
    *
    * @param toPrint
    */
  def toString(toPrint:List[PrintElem]):Unit = {
    val codons:List[(Base, Base, Base)] = getCodons(codonNumb).asInstanceOf[List[Tuple3[Base, Base, Base]]]
    toPrint.foreach(elem => {
      elem match {
        case PrintElem.codeTable => {
          var codeTableAsText = ""
          for (
            codon <- 0 until codons.length;
            aa <- 0 until initAA.length
          ) {
            codeTableAsText += codons(codon).toString() + " -> " + initAA(aa).toString + codeTable(codon)(aa).toString() + "\n"
          }
          "CodeTable\nShows for each codon the possible translations into an amino acid and the aaRSs that enable this.\n" + codeTableAsText
        }
         case PrintElem.mRNA => printMRNA(mRNA, getCodons(codonNumb))
        //case PrintElem.tRNAs => println("tRNAs:\n" + tRNAs.toString()); println()
        //case PrintElem.allTRNA => print2DimTRNAmatrix(allTRNA)
        case PrintElem.`livingAARSs` => livingAARSs.toString()
        case PrintElem.allAARS => print3DimAARSmatrix(allAARS)
      }
    })
  }

  /**
    *
    * @param toHTML
    * @return
    */
  def toHtmlString(toHTML:List[PrintElem]):String = {
    val codons = getCodons(codonNumb)
    var content = "<!DOCTYPE html >\n<html>\n\t<head>\n\t\t<style>\n\t\t\ttable {\n\t\t\t\tfont-family: arial, sans-serif;\n\t\t\t\tborder-collapse: collapse;\n\t\t\t\twidth: 100%;\n\t\t\t}\n\n\t\t\ttd, th {\n\t\t\t\tborder: 1px solid #dddddd;\n\t\t\t\ttext-align: left;\n\t\t\t\tpadding: 8px;\n\t\t\t}\n\n\t\t\ttr:nth-child(even) {\n\t\t\t\tbackground-color: #dddddd;\n\t\t\t}\n\t\t</style>\n\t</head>\n\t<body>"
    content += "</br><h1>Generation " + generationID + "</h1>"
    toHTML.foreach(elem => {
      elem match {
        case PrintElem.codons => {
          content += "<h2>Codons</h2>\n"
          content += codons.mkString("<p>",", ","</p>") + "</br>\n"
        }
        case PrintElem.mRNA => {
          content += "<h2>mRNA</h2>\n"
          content += mRNAtoHTML(mRNA, getCodons(codonNumb))
        }
        case PrintElem.livingAARSs => {
          content += "<h2>Living aaRS</h2>\n<table> \n <tr>\n<th>aaSeq</th>\n<th>Lifeticks</th>\n<th>translations: Amino Acid, Anticodon/Codon, tRNA stem</th>\n</tr>"
          livingAARSs.foreach(aars =>{
            content += "<tr>\n<td>" + aars.aaSeq.mkString(", ") + "</td>\n<td>" + aars.lifeticks + "</td><td>" + aars.translationsToHtmlString() + "</td>\n</tr>"
          })
          content += "</table>"
        }
        case PrintElem.allAARS => print3DimAARSmatrix(allAARS)
        case PrintElem.codeTable => {
          content += "<h2>Code Table</h2>\n<table>\n"
          content += "<tr>\n<th></th>"
          for(
            aa <- 0 until initAA.length
          ){
            content += "<th>" + initAA(aa).toString() + "</th>\n"
          }
          content += "\n"

          for(
            codonPos <- 0 until codonNumb
          ){
            content += "</tr>\n<tr>\n<td>" + codons(codonPos) + "</td>\n"
            for(
              aa <- 0 until initAA.length
            ){
              val aaRSList:List[AARS] = codeTable(codonPos)(aa)
              if(aaRSList != null){
              content += "<td>"
              for(
                aaRS <- aaRSList
              ){
                val r = aaRS.aaSeq(0).id * initAA.length
                val g = aaRS.aaSeq(1).id * initAA.length
                val b = aaRS.aaSeq(2).id * initAA.length
                val rback = 255 - aaRS.aaSeq(0).id * initAA.length
                val gback = 255 - aaRS.aaSeq(1).id * initAA.length
                val bback = 255 - aaRS.aaSeq(2).id * initAA.length
                content += s"<p style='color: rgb($r, $g, $b); background-color: rgb($rback, $gback, $bback);'>" + aaRS.aaSeq.mkString(", ") + "</p></br>\n"
              }
              content += "</td>"
              }else{
                content += "<td></td>"
              }
            }
          }
          content += "\n</tr>\n</table>\n</body>\n</html>"
        }
        //case PrintElem.tRNAs => println("tRNAs:\n" + tRNAs.toString()); println()
        //case PrintElem.allTRNA => print2DimTRNAmatrix(allTRNA)
      }
    })
    content
  }








}