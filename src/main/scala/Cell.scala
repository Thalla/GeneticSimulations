import AA.AA
import Base.Base
import SimulationData.protocol
import Earth._
import PrintElem.PrintElem

import scala.collection.mutable.{ArrayBuffer, ListBuffer}
import scala.specialized
import scala.util.Random


case class CellData(unambiguousness:Double, mRNAdata:ListBuffer[Int], numbAaWithTransl:Int, aaHasTransl:Array[Boolean])
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
class Cell (var livingAARSs:ListBuffer[AARS]) {

  var generationID:Int = 0
  def maxAnticodonNumb:Int = 6
  def r = new scala.util.Random(22)
  var codeTable:Array[Array[ListBuffer[AARS]]] = Array.ofDim[ListBuffer[AARS]](codonNumb,aaNumb)
  //var codeTable2:Array[ListBuffer[AARS]] = new Array[ListBuffer[AARS]](codonNumb*aaNumb)

  /**
    *
    * @return
    */
  def translate():CellData = {

    // data init
    var unambiguousness:Double = 0.0        //Array.fill(codonNumb)(0.0)  //Eindeutigkeit
    var mRNAdata:ListBuffer[Int] = new ListBuffer()

    val (numbAaWithTransl, aaHasTransl) = updateCodeTable()
    this.codeTable = codeTable

    //update lifeticks
    var newLivingAARSs:ListBuffer[AARS] = ListBuffer()
    val livingAARSIt = livingAARSs.iterator
    while (livingAARSIt.hasNext){
      val aaRS = livingAARSIt.next().reduceLifeTicks()
      if (aaRS.lifeticks > 0) newLivingAARSs += aaRS
    }

    val mrnaIt = mRNA.iterator
    while (mrnaIt.hasNext){
      val gene = mrnaIt.next()
      var sequence = new ListBuffer[AA.Value]
      var isStopped: Boolean = false
      val geneIt = gene.iterator
      while(geneIt.hasNext && !isStopped){
        val codon = geneIt.next()
        //var aaRSCounter = 0
        var translationCounter = 0
        var bestAARS: AARS = null
        var max: Double = 0.0
        var bestTranslation: Tuple2[AA, Int] = (null, 0)
        //find aaRS with best translation considering the translation fitness
        for (
          aaRSs: ListBuffer[AARS] <- codeTable(codon) //find all aaRS that can translate this codon
        ) {
          if (!aaRSs.isEmpty) {
            translationCounter += 1
            for (
              aaRS <- aaRSs
            ) {
              //aaRSCounter += 1
              val translationKeys = aaRS.translations.keySet
              //find translation with best fitness of all possible translations
              for (
                tk <- translationKeys
              ) {
                if(tk._2 == codon){
                  val prob: Tuple2[Double, Int] = (aaRS.translations(tk)(0)) //TODO what if more than one item in list?? This is the case when the translation from codon to aa is the same but the codon can have different stems and therefore different probabilities
                  if (max < prob._1) { // what if == ? -> make list and choose randomly from it TODO
                    bestTranslation = tk
                    max = prob._1
                    bestAARS = aaRS
                  }
                }
              }
            }
          }
        }
        if(translationCounter != 0){
          sequence +=  bestTranslation._1
          mRNAdata += (bestTranslation._1.id)

          if (codon < codonNumb) {
            unambiguousness += (1.0/translationCounter.toDouble)/codonNumb.toDouble  //aaRSCounter.toDouble
          }
        }else{
          isStopped = true
          /*var l = sequence.length
            while(l < 3){
              mRNAdata = "NA" :: mRNAdata
              l += 1
            }*/
        }
        //aaRSCounter = 0
      }
      if (!isStopped) {
        val newAARS = allAARS(sequence.remove(0).id)(sequence.remove(0).id)(sequence.remove(0).id)
        newAARS.resetLifeticks()
        if (!newLivingAARSs.contains(newAARS)) {
          newLivingAARSs += newAARS
        }
        sequence = ListBuffer()
      }
    }
    generationID += 1
    this.livingAARSs = newLivingAARSs
    new CellData(unambiguousness, mRNAdata, numbAaWithTransl, aaHasTransl)
  }


  def updateCodeTable():(Int, Array[Boolean]) ={
    val newAaHasTransl = new Array[Boolean](20)
    var translatedAaCounter = 0
    var i = 0
    //initialize
        while (i < aaNumb){
          newAaHasTransl(i) = false
          i += 1
        }


    //new codeTable
    val livingAarsIt = livingAARSs.iterator
    while (livingAarsIt.hasNext){
      val aaRS = livingAarsIt.next()
      val aarsTranslIt = aaRS.translations.iterator
      while (aarsTranslIt.hasNext){
        val translation = aarsTranslIt.next()
        if (newAaHasTransl(translation._1._1.id) == false){
          newAaHasTransl(translation._1._1.id) = true
          translatedAaCounter += 1
        }
        if(codeTable(translation._1._2)(translation._1._1.id) == null) codeTable(translation._1._2)(translation._1._1.id) = ListBuffer(aaRS)
        else codeTable(translation._1._2)(translation._1._1.id) += aaRS
      }
    }
    (translatedAaCounter, newAaHasTransl)
  }

  /**
    *
    * @param sequence
    * @return
    */
  /*def getAARS(sequence:ListBuffer[AA.Value]):AARS={

    var newAARS: AARS = allAARS(sequence.remove(0).id)(sequence.remove(1).id)(sequence.remove(2).id)
    if (newAARS != null) {
      newAARS.resetLifeticks()
    }
    else { //currently not used because all aaRS are created yet
      var translations:Map[(AA, Int),List[(Double, Int)]] = Map()
      val possibleAnticodonNumbs = List.range(1,maxAnticodonNumb)
      val numbAnticodons= possibleAnticodonNumbs(r.nextInt(maxAnticodonNumb))
      for(
        _ <- 0 until numbAnticodons
      ){
        translations += (initAA(r.nextInt(initAA.length)), r.nextInt(codonNumb))-> List((1.0,r.nextInt(codonNumb))) //TODO choose similar codons (first random, rest +x)
      }
      newAARS = new AARS(sequence, translations)
      allAARS(sequence(0).id)(sequence(1).id)(sequence(2).id) = newAARS

    }
    newAARS
  }*/


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

/*
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
  }*/

  /**
    *
    * @param toHTML
    * @return
    */
  def toHtmlString(toHTML:List[PrintElem]):String = {
    val codons = getCodons(codonNumb)
    var content = "<!DOCTYPE html >\n<html>\n\t<head>\n\t\t<style>\n\t\t\ttable {\n\t\t\t\tfont-family: arial, sans-serif;\n\t\t\t\tborder-collapse: collapse;\n\t\t\t\twidth: 100%;\n\t\t\t}\n\n\t\t\ttd, th {\n\t\t\t\tborder: 1px solid #dddddd;\n\t\t\t\ttext-align: left;\n\t\t\t\tpadding: 4px; margin: 0px;\n\t\t\t}\n\n\t\t\ttr:nth-child(even) {\n\t\t\t\tbackground-color: #dddddd;\n\t\t\t}\n\t\t</style>\n\t</head>\n\t<body>"
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
              val aaRSList:ListBuffer[AARS] = codeTable(codonPos)(aa)
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
                content += s"<span style='color: rgb($r, $g, $b); background-color: rgb($rback, $gback, $bback);'>" + aaRS.aaSeq.mkString(", ") + "</span></br>\n"
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