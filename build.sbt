name := "Simple0"

version := "0.1"

scalaVersion := "2.12.6"

libraryDependencies += "org.scala-lang.modules" %% "scala-parser-combinators" % "1.1.0"
//libraryDependencies += "org.biojava" % "biojava3-core" % "3.0"
//libraryDependencies += "org.scala-lang" % "scala-swing" % "2.12"
libraryDependencies += "org.scala-lang.modules" %% "scala-swing" % "2.0.0-M2"
libraryDependencies += "org.ddahl" %% "rscala" % "3.2.2"
unmanagedJars in Compile += file("C:\\Users\\feroc\\OneDrive\\Dokumente\\R\\win-library\\3.5\\rscala\\java\\scala-2.11\\rscala_2.11-3.2.2.jar")

