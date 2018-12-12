import java.io.File

import AA.AA
import PrintElem.PrintElem
import com.github.tototoshi.csv.{CSVReader, CSVWriter}

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
                }
            aa = initAA(i)
          }

        if(translationCounter != 0){
          sequence +=  aa

          if (codon < codonNumb) {
            unambiguousness += (1.0/translationCounter.toDouble)/codonNumb.toDouble  //aaRSCounter.toDouble
          }
        }else{
          isStopped = true
        }
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


  /**
    *
    */
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


  /** creates a start cell, each cell is one generation
    * @param aaRSnumb Number of how many different aaRS proteins exist initially. (see old scaladoc for notes)
    * @param aarsLength The number of amino acids an aaRS consists of. TODO: aaRS can have varying lengths
    * @param initAA initially existing and used amino acids TODO: allow having new/non proteinogenic amino acids TODO start with non essential amino acids
    * @param codonNumb Either 16 (set of twoTuples is created) or 64 (set of real codons is created) or 48 (set of strong commafree codons is created) (this version uses 64, others are not tested)
    * @return start cell
    */
  def init (path:String, aaRSnumb:Int = 22, aarsLength:Int = 3, initAA:Vector[AA] = AA.values.toVector, codonNumb:Int = 64, newLivingAARS:Boolean = false, newAARS:Boolean = false, similarAARS:Boolean = true, newMRNA:Boolean = false):Unit= {
    this.codonNumb = codonNumb
    this.initAA = initAA
    this.aaNumb = initAA.length
    this.path = path
    // create codons
    val codons = getCodons(codonNumb)
    allAARS = new Array[Array[Array[AARS]]](aaNumb)
    codeTable = Array.ofDim[Double](codonNumb, aaNumb)

    // create mRNA (List of genes. Each gene is a List of codon IDs as long as aarsLength)
    if(newMRNA){
      val mRNA:List[List[Int]] = getRandomMRNA(codons.toList, aarsLength, aaRSnumb)
      writeMRNAtoFile(codons.toList, aarsLength, aaRSnumb, new File(path+"mRNA.csv"))
    }

    // create aaRS file
    if(newAARS == true){
      if(similarAARS == true){
        writeAARSwithSimilarityToFile(codons.toList, aarsLength, aaRSnumb, codonNumb, initAA, new File(path+"aaRS.csv"))
      } else {
        writeAARStoFile(codons.toList, aarsLength, aaRSnumb, codonNumb, initAA, new File(path+"aaRS.csv"))
      }
    }

    if(newLivingAARS == true){
      //create living aaRS
      writeLivingAarsToFile(aaRSnumb)
    }

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

  /** One mRNA holds genes for all initially existing aaRS (see aaRSnumb).
    * Gene length is same as aaRS length. No gene exists twice.
    *
    * @param codons codons/twoTuples from which is randomly chosen to build the mRNA
    * @param geneLength number of codons/twoTuples that code for one aaRS
    * @param geneNumb number of aaRS
    * @return mRNA having genes having codons; for example: Vector(Vector(AA, CU, AG), Vector(UA, CC, GC),...)
    */
  def getRandomMRNA(codons:List[Any], geneLength:Int, geneNumb:Int):List[List[Int]] = {
    var mRNA:List[List[Int]] = List()
    val r = new scala.util.Random(22)
    val getRandomGene = () => {
      var gene:List[Int] = List()
      for {
        _ <- 0 to (geneLength-1)
      } gene = gene :+ r.nextInt(codons.length)
      gene
    }
    // get a gene that isn't already existing
    val getUniqueRandomGene = () => {
      var gene:List[Int] = getRandomGene()
      while(mRNA.contains(gene)){
        gene = getRandomGene()
      }
      gene
    }
    //create mRNA with genes for desired start number of aaRS
    for {
      _ <- 0 to (geneNumb-1)
    } mRNA = mRNA :+ getUniqueRandomGene()

    mRNA
  }

  /**
    * not random
    * @param codons
    * @param geneLength
    * @param geneNumb
    * @param file
    */
  def writeMRNAtoFile(codons:List[Any], geneLength:Int, geneNumb:Int, file:File):Unit = {
    var mRNA:List[List[String]] = List()
    val writer = CSVWriter.open(file)
    val r = new scala.util.Random(22)
    var codonCounter = 0
    for(
      _ <- 0 until geneNumb
    ){
      var gene:List[String] = List()
      for {
        _ <- 0 to (geneLength-1)
      } {
        if(codonCounter < codons.length) {
          gene = gene :+ codonCounter.toString()
          codonCounter += 1
        }else{ //each codon was used once
          gene = gene :+ r.nextInt(codons.length).toString()
        }
      }
      mRNA = mRNA :+ gene
    }

    writer.writeAll(mRNA)
    writer.close()
  }

  /**
    *
    * @param codons
    * @param geneLength
    * @param geneNumb
    * @param file
    */
  def writeAARStoFile(codons:List[Any], aarsLength:Int, aaRSnumb:Int, codonNumb: Int, initAA:Vector[AA], file:File):Unit = {
    CSVWriter.open(file).close()
    val writer = CSVWriter.open(file, append= true)
    val r = new scala.util.Random(22)
    // The current genetic code doesn't have more than six different codons per amino acid. The start cell has the same restrictions.
    val numbOfAnticodonsForOneAARS = List(1,2,3,4,5,6)  //TODO remove redundancy: this is the same as maxAnticodonNumb in Cell

    //create 20^3 aaRS
    for(
      i <- 0 until initAA.length;
      j <- 0 until initAA.length;
      k <- 0 until initAA.length
    ){
      val aaRSline:String = ""
      //var translations:Map[(AA, Int),List[(Double, Int)]] = Map() //each aaRS has a translation table that maps amino acids and codons to a probability and tRNA stems. Currently the probability is not used and always 1.0.
      var antiCodonNumb = numbOfAnticodonsForOneAARS(r.nextInt(numbOfAnticodonsForOneAARS.length))      //number of anticodons that are recognized by the aaRS is chosen randomly

      for(
        _ <- 0 until antiCodonNumb
      ){
        // aa1 aa2 aa3, aa, anticodon, probability, stem
        writer.writeRow(Seq(i, j, k, r.nextInt(initAA.length), r.nextInt(codons.length), r.nextDouble().toString(), r.nextInt(codons.length) ))
      }
    }
    writer.close()
  }

  /**
    *
    * @param codons
    * @param geneLength
    * @param geneNumb
    * @param file
    */
  def writeAARSwithSimilarityToFile(codons:List[Any], aarsLength:Int, aaRSnumb:Int, codonNumb: Int, initAA:Vector[AA], file:File):Unit = {
    CSVWriter.open(file).close()
    val writer = CSVWriter.open(file, append= true)
    val r = new scala.util.Random(22)
    // The current genetic code doesn't have more than six different codons per amino acid. The start cell has the same restrictions.
    val numbOfAnticodonsForOneAARS = List(1,2,3,4,5,6)  //TODO remove redundancy: this is the same as maxAnticodonNumb in Cell

    //create 20^3 aaRS
    for(
      i <- 0 until initAA.length
    ){
      val aa = r.nextInt(initAA.length) // every aaRS amino acid sequence that starts with the same amino acid translates to the same amino acid. TODO the same aa as it starts with??
      for (
        j <- 0 until initAA.length

      ){
        val codon = r.nextInt(codonNumb)
        for (
          k <- 0 until initAA.length
        ){
          val aaRSline: String = ""
          //var translations:Map[(AA, Int),List[(Double, Int)]] = Map() //each aaRS has a translation table that maps amino acids and codons to a probability and tRNA stems. Currently the probability is not used and always 1.0.
          var antiCodonNumb = numbOfAnticodonsForOneAARS(r.nextInt(numbOfAnticodonsForOneAARS.length)) //number of anticodons that are recognized by the aaRS is chosen randomly

          for (
            _ <- 0 until antiCodonNumb
          ){
            // aa1 aa2 aa3, aa, anticodon, probability, stem
            val rand = r.nextInt()
            if(rand >= 0.5){
              writer.writeRow(Seq(i, j, k, aa, codon, r.nextDouble().toString(), r.nextInt(codons.length)))
            } else if (rand >= 0.25){
              writer.writeRow(Seq(i, j, k, r.nextInt(initAA.length), codon, r.nextDouble().toString(), r.nextInt(codons.length)))
            } else {
              writer.writeRow(Seq(i, j, k, aa, r.nextInt(codonNumb), r.nextDouble().toString(), r.nextInt(codons.length)))
            }
          }
        }
      }
    }
    writer.close()
  }

  /**
    * create living aaRS file (create random aaSeqences that will come to life)
    * @return
    */
  def writeLivingAarsToFile(aaRSnumb:Int):Unit ={
    val file = new File(path+"livingAARS.csv")
    CSVWriter.open(file).close()
    val writer = CSVWriter.open(file, append = true)
    for(
      _ <- 0 until aaRSnumb
    ){
      writer.writeRow(Seq(r.nextInt(aaNumb), r.nextInt(aaNumb), r.nextInt(aaNumb)))
    }
    writer.close()
  }



}