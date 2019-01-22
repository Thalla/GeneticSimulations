import java.io.File

import AA._
import PrintElem.PrintElem
import com.github.tototoshi.csv.{CSVReader, CSVWriter}

import scala.collection.mutable.ListBuffer
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
class Cell (val r:Random) {
  var path = ""

  var simulationData:SimulationData = null
  var toPrint:List[PrintElem] = List()

  var codonNumb:Int = 0
  var mRNA:List[List[Int]] = List()
  var aaNumb:Int = _
  var initAA:Vector[AA] = Vector()
  private[this] var _aaRSnumb: Int = _
  def aaRSnumb: Int = _aaRSnumb
  def aaRSnumb_=(value: Int): Unit = {
    _aaRSnumb = value
  }
  var allAars:Array[Array[Array[AARS]]] = Array()
  var livingAars:ListBuffer[AARS] = ListBuffer()

  var generationID:Int = -1

  var codeTable = Array.ofDim[Double](codonNumb, aaNumb)
  private[this] var _codeTableFitness: Double = _
  def codeTableFitness: Double = _codeTableFitness
  private[this] def codeTableFitness_=(value: Double): Unit = {
    _codeTableFitness = value
  }

  private[this] var _aaTranslData: (Int, Array[Boolean]) = (0,Array.fill[Boolean](20)(false))  //the second variable exists to remember information from the previous generation to calculate aaChanges
  var aaChanges:Array[Int] = new Array[Int](2)   //genIds, added amino acids, removed amino acids
  def aaTranslData: (Int, Array[Boolean]) = _aaTranslData

  private[this] var _unambiguousness: Double = _
  private[this] def unambiguousness: Double = _unambiguousness
  private[this] def unambiguousness_=(value: Double): Unit = {
    _unambiguousness = value
  }

   var time1:Long = 0
  var time2:Long = 0
  var time0:Long = 0
  var translTime: Long = 0
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
    val livingAARSIt = livingAars.iterator
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

          /*if (codon < codonNumb) {
            unambiguousness += (1.0/translationCounter.toDouble)/codonNumb.toDouble  //aaRSCounter.toDouble
          }*/
        }else{
          isStopped = true
        }
      }
      if (!isStopped) {
        val newAARS = allAars(sequence.remove(0).id)(sequence.remove(0).id)(sequence.remove(0).id)
        newAARS.resetLifeticks()
        if (!newLivingAARSs.contains(newAARS)) {
          newLivingAARSs += newAARS
        }
        sequence = ListBuffer()
      }
    })
    /*codeTableFitness = (aaTranslData._1.toDouble/aaNumb.toDouble)* unambiguousness*/
    this.livingAars = newLivingAARSs
  }


  /**
    *
    */
  def updateCodeTable():Unit ={
    val newAaHasTransl = Array.fill[Boolean](20)(false)
    var translatedAaCounter = 0
    val codonTransl = Array.fill[Boolean](codonNumb*aaNumb)(false)
    val codonTranslCounter = Array.fill[Int](codonNumb)(0)

    //new codeTable
    codeTable = Array.ofDim[Double](codonNumb, aaNumb)
    val livingAarsIt = livingAars.iterator
    while (livingAarsIt.hasNext){
      val aaRS = livingAarsIt.next()
      val aarsTranslIt = aaRS.translations.iterator
      while (aarsTranslIt.hasNext){
        val translation = aarsTranslIt.next()
        if (!newAaHasTransl(translation._1._1.id)){
          newAaHasTransl(translation._1._1.id) = true
          translatedAaCounter += 1
        }
        val fieldPos = translation._1._2*(aaNumb)+translation._1._1.id
        if(!codonTransl(fieldPos)){
          codonTransl(fieldPos) = true
          codonTranslCounter(translation._1._2) += 1
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

    var sum = 0.0
    val cNumb = codonNumb.toDouble
    val lit = codonTranslCounter.iterator
    while (lit.hasNext) {
      val elem = lit.next().toDouble
      val n:Double = 0.0
      if(elem != n){
        sum +=  (1.0/elem)/cNumb
      }

    }
    unambiguousness = sum

    codeTableFitness = (translatedAaCounter.toDouble/aaNumb.toDouble)* unambiguousness


  }
  def initProperties(codonNumb: Int, initAaNumb: Int, mrnaPath: String, aarsPath: String, aarsLifeticksStartValue: Int, livingAarsPath: String, livingAarsDir:String, outputSeed: Int, addToOutput:Boolean): (Boolean, String) = {
    this.codonNumb = codonNumb
    val codons = getCodons(codonNumb)
    this.aaNumb = initAaNumb

    codeTable = Array.ofDim[Double](codonNumb, aaNumb) //init codeTable
    mRNA = readMrna(mrnaPath) //init mRNA
    initAA = AA.values.toVector.take(initAaNumb) //init initAA
    allAars = readAars(aarsPath, initAaNumb, aarsLifeticksStartValue) //init aaRS
    livingAars = readLivingAars(livingAarsPath, allAars, aarsLifeticksStartValue) //init livingAARS

    //init outputFolder
    val outputPath = livingAarsDir + s"\\output_s$outputSeed\\"
    if (!new File(outputPath).exists()) {
      new File(outputPath).mkdirs()
      (true, outputPath)
    } else if (addToOutput) {
      var i = 0 // exists in any case because there are already files in this path therefore a set of living aaRS as well.
      while ((new File(outputPath + s"livingAars${i + 1}.csv").exists())) {
        i += 1
      }
      livingAars = readLivingAars(outputPath + s"livingAars$i.csv", allAars, aarsLifeticksStartValue)
      if (livingAars.length != 0) (true, outputPath)
      else (false, outputPath) // The path already exists and addToOutput is true but there is no file with the last living aaRS set from the last simulation so nothing to build upon.
    } else (false, outputPath) // The path already exists. This means that this param combination has been simulated already. addToOutput is false so no second set of results shall be added to the already existing one.
  }



  /** initializes the cell and the file structure
    * @param geneNumb Number of how many different aaRS proteins exist initially. (see old scaladoc for notes)
    * @param geneLength The number of amino acids an aaRS consists of. TODO: aaRS can have varying lengths
    * @param initAA initially existing and used amino acids TODO: allow having new/non proteinogenic amino acids TODO start with non essential amino acids
    * @param codonNumb Either 16 (set of twoTuples is created) or 64 (set of real codons is created) or 48 (set of strong commafree codons is created) (this version uses 64, others are not tested)
    */
  def init (basePath:String, mrnaSeed:Int, codonNumb:Int, geneLength:Int, geneNumb:Int, mrnaId:Int, initAaNumb:Int, similarAars:Boolean = true, aarsSeed:Int, maxAnticodonNumb:Int, aarsLifeticksStartValue:Int, livingAarsStartNumb:Int, livingAarsSeed:Int,  outputSeed:Int, addToOutput:Boolean):Boolean= {
    path = basePath

    this.codonNumb = codonNumb
    val codons = getCodons(codonNumb)
    this.aaNumb = initAaNumb

    //init codeTable
    codeTable = Array.ofDim[Double](codonNumb, aaNumb)

    //init mRNA
      val mrnaName = writeNewMrna(path, mrnaSeed, codons, geneLength, geneNumb, mrnaId)
      path += s"$mrnaName\\"
      mRNA = readMrna(path + s"$mrnaName.csv")

    //init initAA
    initAA = AA.values.toVector.take(initAaNumb)
    path += s"initAaNumb_$initAaNumb\\"

    //init aaRS
    path += s"similar_$similarAars\\"
    val aarsName = writeNewAars(path, similarAars, aarsSeed, codonNumb, initAaNumb, maxAnticodonNumb, aarsLifeticksStartValue)
    path += s"$aarsName\\"
    allAars = readAars(path + s"$aarsName.csv", initAaNumb, aarsLifeticksStartValue)

    //init livingAARS
    val livingAarsName = writeNewLivingAars(path, livingAarsSeed, livingAarsStartNumb, initAaNumb)
    path += s"$livingAarsName\\"
    livingAars = readLivingAars(path + s"$livingAarsName.csv", allAars, aarsLifeticksStartValue)

    //init outputFolder
    path += s"output_s$outputSeed\\"
    if(!new File(path).exists()){
      new File(path).mkdirs()
      true
    } else if(addToOutput) {
      var i = 0   // exists in any case because there are already files in this path therefore a set of living aaRS as well.
      while ((new File(path + s"livingAars${i+1}.csv").exists())) {
        i += 1
      }
      livingAars = readLivingAars(path + s"livingAars$i.csv", allAars, aarsLifeticksStartValue)
      if(livingAars.length != 0) true
      else false
    }else false
  }

  /**
    * Creates new mRNA if a mRNA with the given parameters doesn't exist yet.
    * @param basePath
    * @param mrnaSeed
    * @param codonNumb
    * @param initAaNumb
    * @param id
    * @return fileName: name of the mRNA file that fits to the given parameters
    */
  def writeNewMrna(path:String, mrnaSeed:Int, codons:IndexedSeq[(Any)], geneLength:Int, geneNumb:Int, id:Int): (String, Boolean) = {
    val mrnaName = s"mRNA_s$mrnaSeed" + s"_c${codons.length}" + s"_gl$geneLength" + s"_gn$geneNumb" + s"_#$id"
    val mrnaPath = path + mrnaName + "\\"
    var newCombi = false
    if (!new File(mrnaPath).exists()) {
      newCombi = true
      new File(mrnaPath).mkdirs()
      val mRNA: List[List[Int]] = getRandomMRNA(codons.toList, geneLength, geneNumb)
      writeMRNAtoFile(codons.toList, geneLength, geneNumb, new File(mrnaPath + mrnaName + ".csv"))
    }
    (mrnaName, newCombi)
  }

  def readMrna(path:String): List[List[Int]] ={
    var reader = CSVReader.open(new File(path))
    val data = reader.all()
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
    mRNA
  }

  /**
    *
    * @param path
    * @param similarAars
    * @param aarsSeed
    * @param codons
    * @param geneLength
    * @param geneNumb
    * @return
    */
  def writeNewAars(path:String, similarAars:Boolean, aarsSeed:Int, codonNumb:Int, initAaNumb:Int, maxAnticodonNumb:Int, lifeticksStartValue:Int):(String, Boolean) ={
    val aarsName = s"aaRS_s${aarsSeed}_ac${maxAnticodonNumb}_lt$lifeticksStartValue"
    val filePath = path + s"$aarsName\\"
    var newCombi = false
    if( ! new File(filePath).exists) {
      newCombi = true
      new File(filePath).mkdirs()
      if(similarAars) writeAARSwithSimilarityToFile(codonNumb, initAaNumb, maxAnticodonNumb, new File(filePath+ s"$aarsName.csv"), aarsSeed)
      else writeAARStoFile(codonNumb, initAaNumb, maxAnticodonNumb, new File(filePath+ s"$aarsName.csv"), aarsSeed)
    }
    (aarsName, newCombi)
  }

  /**
    *
    * @param path
    * @param initAaNumb
    * @param aarsLifeticksStartValue
    * @return
    */
  def readAars(path:String, initAaNumb:Int, aarsLifeticksStartValue:Int): Array[Array[Array[AARS]]] ={
    allAars = new Array[Array[Array[AARS]]](initAaNumb)
    val reader = CSVReader.open(new File(path))
    val data = reader.all()
    var translations:Map[(AA, Int),List[(Double, Int)]] = Map()
    // Array with all existing aaRS, initialised with placeholders
    //var allAARS:Array[Array[Array[AARS]]] = Array.fill[Array[Array[AARS]]](aaNumb)(Array.fill[Array[AARS]](aaNumb)(Array.fill[AARS](aaNumb)(new AARS(Vector(), translations))))
    // give each aaRS in allAARS the correct aaSeq (is dependent from its position in allAARS)
    for(
      i <- 0 until initAaNumb
    ){
      allAars(i) = new Array [Array[AARS]](initAaNumb)
      for(
        j <- 0 until initAaNumb
      ){
        allAars(i)(j) = new Array [AARS](initAaNumb)
        for(
          k <- 0 until initAaNumb
        ){
          val aaSeq:Vector[AA]= Vector(initAA(i), initAA(j), initAA(k))
          allAars(i)(j)(k) = new AARS(aaSeq, translations, aarsLifeticksStartValue)
        }
      }
    }

    // read file data and give each aaRS its translations
    for(
      line <- data
    ){
      //get Array Index of aaRS
      val i1 = line.head.toInt
      val i2 = line(1).toInt
      val i3 = line(2).toInt
      // add translation to aaRS
      allAars(i1)(i2)(i3).translations += ((initAA(line(3).toInt), line(4).toInt)->List((line(5).toDouble, line(6).toInt)))
    }
    reader.close()

    allAars
  }

  /**
    *
    * @param path
    * @param livingAarsId
    * @param geneNumb
    * @return fileName
    */
  def writeNewLivingAars(path: String, livingAarsSeed: Int, livingAarsNumb: Int, initAaNumb:Int): (String, Boolean) ={
    val livingAarsName = s"livingAars_n${livingAarsNumb}_s$livingAarsSeed"
    val filePath = path + s"$livingAarsName\\"
    var newCombi = false
    if( ! new File(filePath).exists) {
      new File(filePath).mkdirs()
      writeLivingAarsToFile(filePath + s"$livingAarsName.csv", livingAarsSeed, livingAarsNumb, initAaNumb)
    }
    (livingAarsName, newCombi)
  }

  //Mostly same as readLivingAars!!!!!!!!! REDUNDANT
  def readLivingAarsWithLt(path: String, allAARS: Array[Array[Array[AARS]]]): ListBuffer[AARS] ={
    // read the living aaRS
    val f = new File(path)
    if(f.length != 0){
      val reader = CSVReader.open(f)
      val data = reader.all()
      for(
        line <- data
      ){
        allAARS(line(0).toInt)(line(1).toInt)(line(2).toInt).lifeticks = line(3).toInt  //only changed part
        livingAars = livingAars :+ allAARS(line(0).toInt)(line(1).toInt)(line(2).toInt)
      }
      reader.close()
      livingAars = livingAars.distinct
      livingAars
    }
    else ListBuffer()
  }

  def readLivingAars(path: String, allAARS: Array[Array[Array[AARS]]], aarsLifeticksStartValue: Int): ListBuffer[AARS] ={
    // read the living aaRS
    val f = new File(path)
    if(f.length != 0){
      val reader = CSVReader.open(f)
      val data = reader.all()
      for(
        line <- data
      ){
        allAARS(line(0).toInt)(line(1).toInt)(line(2).toInt).lifeticks = aarsLifeticksStartValue
        livingAars = livingAars :+ allAARS(line(0).toInt)(line(1).toInt)(line(2).toInt)
      }
      reader.close()
      livingAars = livingAars.distinct
      livingAars
    }
    else ListBuffer()
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
      case _ => 0 to 3
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
        case 64 | 48 | _ => (Base(i),Base(j), Base(k))
      }})
    codons.take(codonNumb)
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

  def writeRealAARS(codons:List[Any], aarsLength:Int, aaRSnumb:Int, codonNumb: Int, initAA:Vector[AA], file:File, seed:Int):Unit ={
    CSVWriter.open(file).close()
    val writer = CSVWriter.open(file, append= true)
    val numbOfAnticodonsForOneAARS = List(1,2,3,4,5,6)  //TODO remove redundancy: this is the same as maxAnticodonNumb in Cell
    val r = new scala.util.Random(seed)
  }


  /**
    *
    * @param codons
    * @param geneLength
    * @param geneNumb
    * @param file
    */
  def writeAARStoFile(codonNumb:Int, initAaNumb:Int, maxAnticodonNumb:Int, file:File, seed:Int):Unit = {
    CSVWriter.open(file).close()
    val writer = CSVWriter.open(file, append= true)
    //val r = new scala.util.Random(22)
    // The current genetic code doesn't have more than six different codons per amino acid. The start cell has the same restrictions.
    val r = new scala.util.Random(seed)
    //create 20^3 aaRS
    for(
      i <- 0 until initAaNumb;
      j <- 0 until initAaNumb;
      k <- 0 until initAaNumb
    ){
      //var translations:Map[(AA, Int),List[(Double, Int)]] = Map() //each aaRS has a translation table that maps amino acids and codons to a probability and tRNA stems. Currently the probability is not used and always 1.0.
      var antiCodonNumb = r.nextInt(maxAnticodonNumb)+1      //number of anticodons that are recognized by the aaRS is chosen randomly

      for(
        _ <- 0 until antiCodonNumb
      ){
        // aa1 aa2 aa3, aa, anticodon, probability, stem
        writer.writeRow(Seq(i, j, k, r.nextInt(initAaNumb), r.nextInt(codonNumb), r.nextDouble().toString(), r.nextInt(codonNumb) ))
      }
    }
    writer.close()
  }

  /**
    * line from aaRS file:
    * aaSeq ids, aa id, codon id, affinity, ?
    * @param codons
    * @param geneLength
    * @param geneNumb
    * @param file
    */
  def writeAARSwithSimilarityToFile(codonNumb:Int, initAaNumb:Int, maxAnticodonNumb:Int, file:File, seed:Int):Unit = {
    CSVWriter.open(file).close()
    val writer = CSVWriter.open(file, append= true)
    //val r = new scala.util.Random(22)
    // The current genetic code doesn't have more than six different codons per amino acid. The start cell has the same restrictions.
    val r = new scala.util.Random(seed)
    //create 20^3 aaRS
    for(
      i <- 0 until initAaNumb
    ){
      val aa = i//r.nextInt(initAA.length) // every aaRS amino acid sequence that starts with the same amino acid translates to the same amino acid. TODO the same aa as it starts with??
      for (
        j <- 0 until initAaNumb
      ){
        for (
          k <- 0 until initAaNumb
        ){
          val codon = r.nextInt(codonNumb)
          //var translations:Map[(AA, Int),List[(Double, Int)]] = Map() //each aaRS has a translation table that maps amino acids and codons to a probability and tRNA stems. Currently the probability is not used and always 1.0.
          val antiCodonNumb = r.nextInt(maxAnticodonNumb)+1 //number of anticodons that are recognized by the aaRS is chosen randomly

          for (
            _ <- 0 until antiCodonNumb
          ){
            // aa1 aa2 aa3, aa, anticodon, probability, stem
            val rand = r.nextInt()
            if(rand >= 0.8){
              writer.writeRow(Seq(i, j, k, aa, codon, r.nextDouble().toString(), r.nextInt(codonNumb)))
            } else if (rand >= 0.5){
              writer.writeRow(Seq(i, j, k, r.nextInt(initAaNumb), codon, r.nextDouble().toString(), r.nextInt(codonNumb)))
            } else if (rand >= 0.2) {
              writer.writeRow(Seq(i, j, k, aa, r.nextInt(codonNumb), r.nextDouble().toString(), r.nextInt(codonNumb)))
            }
            else {
              writer.writeRow(Seq(i, j, k, r.nextInt(initAaNumb), r.nextInt(codonNumb), r.nextDouble().toString(), r.nextInt(codonNumb)))
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
  def writeLivingAarsToFile(path:String, seed:Int, aaRSnumb:Int, InitAaNumb:Int):Unit ={
    val r = new scala.util.Random(seed)
    val file = new File(path)//+s"livingAARS_#$id.csv")
    CSVWriter.open(file).close()
    val writer = CSVWriter.open(file, append = true)
    for(
      _ <- 0 until aaRSnumb
    ){
      writer.writeRow(Seq(r.nextInt(InitAaNumb), r.nextInt(InitAaNumb), r.nextInt(InitAaNumb)))
    }
    writer.close()
  }



 /* /** creates a start cell, each cell is one generation
    * @param geneNumb Number of how many different aaRS proteins exist initially. (see old scaladoc for notes)
    * @param geneLength The number of amino acids an aaRS consists of. TODO: aaRS can have varying lengths
    * @param initAA initially existing and used amino acids TODO: allow having new/non proteinogenic amino acids TODO start with non essential amino acids
    * @param codonNumb Either 16 (set of twoTuples is created) or 64 (set of real codons is created) or 48 (set of strong commafree codons is created) (this version uses 64, others are not tested)
    * @return start cell
    */
  def init0 (startPath:String, geneNumb:Int = 22, geneLength:Int = 3, aarsLifeticksStartValue:Int = 10, initAaNumb:Int = 20, aarsSeed:Int = 0, maxAnticodonNumb:Int, codonNumb:Int = 64, newLivingAARS:Boolean = false, livingAarsId:Int, newAARS:Boolean = false, similarAARS:Boolean = true, newMRNA:Boolean = false):Unit= {
    this.path = startPath

    initAA = AA.values.toVector.take(initAaNumb)
    allAars = new Array[Array[Array[AARS]]](initAaNumb)
    this.codonNumb = codonNumb
    val codons = getCodons(codonNumb)   // create codons
    codeTable = Array.ofDim[Double](codonNumb, initAaNumb)

    // create mRNA (List of genes. Each gene is a List of codon IDs as long as aarsLength)
    if(newMRNA){
      val mRNA:List[List[Int]] = getRandomMRNA(codons.toList, geneLength, geneNumb)
      writeMRNAtoFile(codons.toList, geneLength, geneNumb, new File(path+"mRNA.csv"))
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

    path += s"$codonNumb\\$initAaNumb\\"

    // create aaRS file
    if(newAARS == true){
      if(similarAARS == true){
        path += s"1\\aaRS_$aarsSeed\\"
        new File(path).mkdirs()
        writeAARSwithSimilarityToFile(codons.toList, geneLength, geneNumb, codonNumb, initAA, maxAnticodonNumb, new File(path+s"aaRS_s$aarsSeed.csv"), aarsSeed)
      } else {
        path += s"0\\aaRS_s$aarsSeed\\"
        new File(path).mkdirs()
        writeAARStoFile(codons.toList, geneLength, geneNumb, codonNumb, initAA, maxAnticodonNumb, new File(path+s"aaRS_s$aarsSeed.csv"), aarsSeed)
      }
    }

    // read aaRS
    reader = CSVReader.open(new File(path+s"aaRS_s$aarsSeed.csv"))
    data = reader.all()
    var translations:Map[(AA, Int),List[(Double, Int)]] = Map()
    // Array with all existing aaRS, initialised with placeholders
    //var allAARS:Array[Array[Array[AARS]]] = Array.fill[Array[Array[AARS]]](aaNumb)(Array.fill[Array[AARS]](aaNumb)(Array.fill[AARS](aaNumb)(new AARS(Vector(), translations))))
    // give each aaRS in allAARS the correct aaSeq (is dependent from its position in allAARS)
    for(
      i <- 0 until initAaNumb
    ){
      allAars(i) = new Array [Array[AARS]](initAaNumb)
      for(
        j <- 0 until initAaNumb
      ){
        allAars(i)(j) = new Array [AARS](initAaNumb)
        for(
          k <- 0 until initAaNumb
        ){
          val aaSeq:Vector[AA]= Vector(initAA(i), initAA(j), initAA(k))
          allAars(i)(j)(k) = new AARS(aaSeq, translations, aarsLifeticksStartValue)
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
      allAars(i1)(i2)(i3).translations += ((initAA(line(3).toInt), line(4).toInt)->List((line(5).toDouble, line(6).toInt)))
    }
    reader.close()




    if(newLivingAARS == true){
      path += s"\\livingAars_#$livingAarsId\\"
      new File(path).mkdirs()
      //create living aaRS
      writeLivingAarsToFile(path, livingAarsId, geneNumb, aaNumb)
    }

    // read the living aaRS
    reader = CSVReader.open(new File(s"livingAars_#$livingAarsId.csv"))
    data = reader.all()
    for(
      line <- data
    ){
      allAars(line(0).toInt)(line(1).toInt)(line(2).toInt).lifeticks = aarsLifeticksStartValue
      livingAars = livingAars :+ allAars(line(0).toInt)(line(1).toInt)(line(2).toInt)
    }
    reader.close()
    livingAars = livingAars.distinct
    //new Cell(mRNA, livingAARSs, allAARS, initAA, codonNumb,0)

  }
*/
}