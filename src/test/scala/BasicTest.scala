import org.junit._
import org.scalatest.junit.AssertionsForJUnit

class BasicTest extends AssertionsForJUnit {

  val path = "C:\\Users\\feroc\\OneDrive\\Dokumente\\HS\\Semester\\4\\Thesis\\Modeling\\testCSVs\\"

  val cell = new Cell()
  @Before def init():Unit={
    cell.init(path, 5, 3, Vector(AA.Phe, AA.Leu), 16, false, false, false)
  // Problem: ich verwende die IDs der Aminosäuren, Ala hat ID 16, aber im allAARS Array hat Ala Position 0, wenn Ala an nullter Stelle im initAA Vektor steht.
  }
  @Test def basics ():Unit = {
    assert(cell.aaNumb == 2)
    cell.translate()
    assert(cell.allAARS(1)(0)(0).lifeticks == 10)
    assert(cell.allAARS(0)(1)(0).lifeticks == 9)
    assert(cell.codeTableFitness == 0.6875)
    cell.allAARS.foreach(x => x.foreach(y => y.foreach(z => println(z.toString()))))
    println(cell.mRNA.toString())
    cell.livingAARSs.foreach(x => println(x.toString()))
  }

}
