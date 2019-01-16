def init0(basePath: (Int, Int) => Int):Int =>  (Int => Int) ={
  def init1(mrnaSeed:Int): Int=>Int ={
    (some:Int) => basePath(mrnaSeed, some)
  }
  init1
}

def x (i: Int, s: Int): Int ={
  8
}

init0(x)

