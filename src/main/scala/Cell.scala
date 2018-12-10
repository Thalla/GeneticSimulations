import java.io.File

import AA.AA
import PrintElem.PrintElem
import com.github.tototoshi.csv.CSVReader

import scala.collection.mutable.{ArrayBuffer, ListBuffer}
import scala.specialized
import scala.util.Random

/** Cell holds the state of the simulation.
  * State: Everything fits to the generation ID except the livingAARS which are the ones of the next generation.
  *        The generation ID changes with the translation. So each generation is one translated set. Except generation -1 .
  * Cell is responsible for computing the next simulation state by translating mRNA and getting a new set of living aaRS
  * TODO connect states, reduce constructor parameters
  * @param mRNA to be translated to get the aaRS and therefore the codeTable for the next cell generation
  * @param codeTable either check the code table or the current aaRSs themselves to find out about the translation from codon to amino acid
  * @param livingAARSs current set of aaRS, used for translating the mRNA to get the aaRS set for the next cell generation
  * @param allAARS this list holds all valid aaRS that could exist and their randomly predefined behaviour
  * @param generationID nextCell has generationID+1
  */
class Cell () {
  var livingAARSs:ListBuffer[AARS] = ListBuffer()
  var generationID:Int = -1
  private[this] var _codeTableFitness: Double = _

  def codeTableFitness: Double = _codeTableFitness

  private[this] def codeTableFitness_=(value: Double): Unit = {
    _codeTableFitness = value
  }
  private [this] val maxAnticodonNumb:Int = 6
  def r = new scala.util.Random(22)
  var codeTable = Array.ofDim[Double](codonNumb, aaNumb)
  var mRNA:List[List[Int]] = List()
  var allAARS:Array[Array[Array[AARS]]] = Array()
  var initAA:Vector[AA] = Vector()
  var codonNumb:Int = 0
  var path = ""
  var aaNumb:Int = _
  private[this] var _aaRSnumb: Int = _

  def aaRSnumb: Int = _aaRSnumb

  def aaRSnumb_=(value: Int): Unit = {
    _aaRSnumb = value
  }

  private[this] var _aaTranslData: (Int, Array[Boolean]) = (0,Array.fill[Boolean](20)(false))  //the second variable exists to remember information from the previous generation to calculate aaChanges
  var aaChanges:Array[Int] = new Array[Int](2)   //genIds, added amino acids, removed amino acids


  def aaTranslData: (Int, Array[Boolean]) = _aaTranslData

  private[this] var _unambiguousness: Double = _

  private[this] def unambiguousness: Double = _unambiguousness

  private[this] def unambiguousness_=(value: Double): Unit = {
    _unambiguousness = value
  }


  //var codeTable:Array[Array[List[AARS]]] = Array.ofDim[List[AARS]](codonNumb,aaNumb)
  //var codeTable2:Array[ListBuffer[AARS]] = Array.fill[ListBuffer[AARS]](codonNumb){new ListBuffer[AARS]()}
  //var codeTable:Array[Array[Int]] = Array.fill[Array[Int]](codonNumb){Array.fill[Int](aaNumb){0}}
   var time1 = 0:Long
  var time2 = 0:Long
  var time0 = 0:Long
  var translTime = 0:Long
  def translTime[R](block: => R): R = {
    val t0 = System.nanoTime()
    val result = block // call-by-name
    val t1 = System.nanoTime()
    //println("Elapsed time: " + (t1 - t0) + "ns")
    translTime += (t1 - t0)
    result
  }
  def time0[R](block: => R): R = {
    val t0 = System.nanoTime()
    val result = block // call-by-name
    val t1 = System.nanoTime()
    //println("Elapsed time: " + (t1 - t0) + "ns")
    time0 += (t1 - t0)
    result
  }
  def time1[R](block: => R): R = {
    val t0 = System.nanoTime()
    val result = block // call-by-name
    val t1 = System.nanoTime()
    //println("Elapsed time: " + (t1 - t0) + "ns")
    time1 += (t1 - t0)
    result
  }
  def time2[R](block: => R): R = {
    val t0 = System.nanoTime()
    val result = block // call-by-name
    val t1 = System.nanoTime()
    //println("Elapsed time: " + (t1 - t0) + "ns")
    time2 += (t1 - t0)
    result
  }


  /**
    *
    * @return
    */
  def translate():Unit= {
    generationID += 1
    unambiguousness = 0.0

    // data init
    var mRNAdata:ListBuffer[Int] = new ListBuffer()

    def timer[R](block: => R): R = {
      val t0 = System.nanoTime()
      val result = block // call-by-name
      val t1 = System.nanoTime()
      println("Elapsed time: " + (t1 - t0) + "ns")
      result
    }


   updateCodeTable()
     //timer(for(_ <- 0 until 5000000) y())
    //timer(for(_<-0 until 5000000)updateCodeTable())

    //update lifeticks
   var newLivingAARSs:ListBuffer[AARS] = ListBuffer()
    val livingAARSIt = livingAARSs.iterator
    //1s
    while (livingAARSIt.hasNext){
      val aaRS = livingAARSIt.next().reduceLifeTicks()
      if (aaRS.lifeticks > 0) newLivingAARSs += aaRS
    }


    val mrnaIt = mRNA.iterator
    // foreach gene
    translTime(  while (mrnaIt.hasNext){
      val gene = mrnaIt.next()
      //start new sequence
      var sequence = new ListBuffer[AA.Value]
      var isStopped: Boolean = false
      val geneIt = gene.iterator
      var genePos = 0 //counts genes
      //foreach codon in gene
      while(geneIt.hasNext && !isStopped){
        val codon = geneIt.next()
        genePos += 1
        //var aaRSCounter = 0
        var translationCounter = 0
        var aa:AA=null
        //find aaRS with best translation considering the translation fitness
        val translations = codeTable(codon)
          if(translations != null){
              //find translation with best fitness of all possible translations
            /*var sum = 0.0
            var counter = 0.0
            val lit = translations.iterator
            var maxID = 0
            var maxValue = 0.0
            var i = 0
            while (lit.hasNext) {
              val value = lit.next()
              sum += value
              if(value > maxValue){
                maxID = i
                maxValue = value
              }
              if(translations(i) != null){
                counter += translations(i)
              }
              i+=1
            }
            while (lit.hasNext) {
              val value = lit.next()
              sum += value
              if(value > maxValue){
                maxID = i
                maxValue = value
              }
              if(translations(i) != null){
                counter += translations(i)
              }
              i+=1
            }
            aa = initAA(maxID)
            val random = r.nextDouble()*13*/

                var sum = 0.0
                val lit = translations.iterator
                while (lit.hasNext) {
                  sum += lit.next()
                }
                val random = r.nextDouble()*sum //random number between 0 and sum
                var counter = 0.0
                var i = -1
                while(counter <= random && i < aaNumb-1){
                  i += 1
                 val n :Double = 0.0
                  if(translations(i) != n){
                    counter += translations(i)
                    translationCounter += 1
                  }
                                  //counter +=  translations.lift(i).getOrElse(0.0)
                }
            aa = initAA(i)
          }

        if(translationCounter != 0){
          sequence +=  aa
          //mRNAdata += (aa.id)

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
    })

    codeTableFitness = (aaTranslData._1.toDouble/aaNumb.toDouble)* unambiguousness
    this.livingAARSs = newLivingAARSs
  }



  def updateCodeTable() ={
    val newAaHasTransl = new Array[Boolean](20)
    var translatedAaCounter = 0
    var i = 0
    //initialize
        while (i < aaNumb){
          newAaHasTransl(i) = false
          i += 1
        }

    //new codeTable
    codeTable = Array.ofDim[Double](codonNumb, aaNumb)
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
        if (codeTable(translation._1._2) != null){
          codeTable(translation._1._2)(translation._1._1.id) = translation._2(0)._1 + codeTable(translation._1._2)(translation._1._1.id)
        } else {
          codeTable(translation._1._2)(translation._1._1.id) = translation._2(0)._1
        }
      }
    }


    aaChanges = new Array[Int](2)
    for(
      i <- 0 until aaNumb
    ){
      (_aaTranslData._2(i), newAaHasTransl(i)) match {
        case (false, true) => aaChanges(0) += 1
        case (true, false) => aaChanges(1) += 1
        case _ =>
      }
    }
    _aaTranslData = (translatedAaCounter, newAaHasTransl)

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

 /*
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

  }*/

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
  /*def toHtmlString(toHTML:List[PrintElem]):String = {
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
              val aaRSList:ListBuffer[AARS] = for(x <- codeTable(codonPos))yield(x._1)
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
  }*/





  /** creates a start cell, each cell is one generation
    * @param aaRSnumb Number of how many different aaRS proteins exist initially. (see old scaladoc for notes)
    * @param aarsLength The number of amino acids an aaRS consists of. TODO: aaRS can have varying lengths
    * @param initAA initially existing and used amino acids TODO: allow having new/non proteinogenic amino acids TODO start with non essential amino acids
    * @param codonNumb Either 16 (set of twoTuples is created) or 64 (set of real codons is created) or 48 (set of strong commafree codons is created) (this version uses 64, others are not tested)
    * @return start cell
    */
  def init (path:String, aaRSnumb:Int = 22, aarsLength:Int = 3, initAA:Vector[AA] = AA.values.toVector, codonNumb:Int = 64):Unit= {
    this.codonNumb = codonNumb
    this.initAA = initAA
    this.aaNumb = initAA.length
    this.aaRSnumb = aaRSnumb
    this.path = path
    // create codons
    val codons = getCodons(codonNumb)
    allAARS = new Array[Array[Array[AARS]]](aaNumb)
    codeTable = Array.ofDim[Double](codonNumb, aaNumb)


    // two mRNA creation versions, first random with no gene twice, second saved in file
    // create mRNA (List of genes. Each gene is a List of codon IDs as long as aarsLength)
    // val mRNA:List[List[Int]] = getRandomMRNA(codons.toList, aarsLength, aaRSnumb)
    // writeMRNAtoFile(codons.toList, aarsLength, aaRSnumb, new File(path+"mRNA.csv"))

    // create aaRS file
    // writeAARStoFile(codons.toList, aarsLength, aaRSnumb, codonNumb, initAA, new File(path+"aaRS.csv"))

    //create living aaRS
    //writeLivingAarsToFile()



    // read mRNA
    var reader = CSVReader.open(new File(path+"mRNA.csv"))
    var data = reader.all()
    //var mRNA:List[List[Int]] = List()
    for(
      line <- data
    ){
      var gene:List[Int] = List()
      for(
        elem <- line
      ){
        gene = gene :+ elem.toInt
      }
      mRNA = mRNA :+ gene
    }
    reader.close()


    // read aaRS
    reader = CSVReader.open(new File(path+"aaRS.csv"))
    data = reader.all()
    var translations:Map[(AA, Int),List[(Double, Int)]] = Map()
    // Array with all existing aaRS, initialised with placeholders
    //var allAARS:Array[Array[Array[AARS]]] = Array.fill[Array[Array[AARS]]](aaNumb)(Array.fill[Array[AARS]](aaNumb)(Array.fill[AARS](aaNumb)(new AARS(Vector(), translations))))
    // give each aaRS in allAARS the correct aaSeq (is dependent from its position in allAARS)
    for(
      i <- 0 until aaNumb
    ){
      allAARS(i) = new Array [Array[AARS]](aaNumb)
      for(
        j <- 0 until aaNumb
      ){
        allAARS(i)(j) = new Array [AARS](aaNumb)
        for(
          k <- 0 until aaNumb
        ){
          val aaSeq:Vector[AA]= Vector(initAA(i), initAA(j), initAA(k))
          allAARS(i)(j)(k) = new AARS(aaSeq, translations)

        }
      }
    }

    // read file data and give each aaRS its translations
    for(
      line <- data
    ){
      //get Array Index of aaRS
      val i1 = line(0).toInt
      val i2 = line(1).toInt
      val i3 = line(2).toInt
      // add translation to aaRS
      allAARS(i1)(i2)(i3).translations += ((initAA(line(3).toInt), line(4).toInt)->List((line(5).toDouble, line(6).toInt)))
    }
    reader.close()


    // read the living aaRS
    val lifeticksStartValue = allAARS(0)(0)(0).lifeticksStartValue
    reader = CSVReader.open(new File(path+"livingAARS.csv"))
    data = reader.all()
    for(
      line <- data
    ){
      allAARS(line(0).toInt)(line(1).toInt)(line(2).toInt).lifeticks = lifeticksStartValue
      livingAARSs = livingAARSs :+ allAARS(line(0).toInt)(line(1).toInt)(line(2).toInt)
    }
    reader.close()

    //new Cell(mRNA, livingAARSs, allAARS, initAA, codonNumb,0)

  }

  /** generate from nucleobases
    * @param codonNumb 16 -> twoTuples, 64 -> codons, 48 -> strong comma free codons
    * @return either list of twoTuples or codons or strong comma free codons
    */
  def getCodons(codonNumb:Int):IndexedSeq[(Any)] = {
    //range for adding third base
    val range:Range = codonNumb  match{
      case 16 => 0 to 0       //-> no third base
      case 64 | 48 => 0 to 3
    }
    // if strong comma free code shall be generated this filter is used
    val filterScf = (k:Int, i:Int) => {
      if(codonNumb == 48) k!=i
      else true
    }

    var codons = (for {
      i <- 0 to 3 // -> four nucleobases
      j <- 0 to 3
      k <- range
      if filterScf(k, i)  // if scf option 48 then filter codons
    } yield {
      codonNumb match{
        case 16 => (Base(i),Base(j))
        case 64 | 48=> (Base(i),Base(j), Base(k))
      }})
    codons
  }

  //print aaRS-Matrix
  def print3DimAARSmatrix (matrix:Array[Array[Array[AARS]]]):Unit = {
    for (
      i <- 0 until matrix.length;
      j <- 0 until matrix(0).length;
      k <- 0 until matrix(0)(0).length
    ) {
      val elem = matrix(i)(j)(k)
      if (elem == null) {}//println(i + ", " + j + ", " + k + " ")
      else println(i + " , " + j + " , " + k + ": " + matrix(i)(j)(k).toString())
    }
  }

  def mRNAtoHTML(mRNA:List[List[Int]], codons:IndexedSeq[(Any)]):String = {
    var content = "<ul>\n"
    mRNA.foreach(gene => {
      content += "<li>\n"
      gene.foreach(codonID => content += codons(codonID).toString() )
      content += "</li>\n"
    })
    content += "</ul>\n"
    content
  }


}