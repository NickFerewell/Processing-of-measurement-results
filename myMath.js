
function tdistr(n, p) { //df, n - degrees of freedom, number of measurements, sl, p - significance level, желаемая точность(probability - вероятность, что ответ правильный, в нужных пределах).
    if(n <= 0 || Math.abs(n) - Math.abs(n.toFixed(0)) != 0 || p <= 0 || p >= 1){
        console.error("Invalid n or p.")
        return 0;
    }
    return subt(n, p);
}

function subt(n, p) {
    if(p == 0.5){
        return 0;
    } else if(p < 0.5){
        return - subt(n, 1 - p);
    }

    var u = _subu(p);
    var u2 = Math.pow(u, 2);

    var a = (u2 + 1) / 4;
            var b = ((5 * u2 + 16) * u2 + 3) / 96;
            var c = (((3 * u2 + 19) * u2 + 17) * u2 - 15) / 384;
            var d = ((((79 * u2 + 776) * u2 + 1482) * u2 - 1920) * u2 - 945)
                / 92160;
            var e = (((((27 * u2 + 339) * u2 + 930) * u2 - 1782) * u2 - 765) * u2
                + 17955) / 368640;

            var x = u * (1 + (a + (b + (c + (d + e / n) / n) / n) / n) / n);

            if (n <= Math.pow(Math.log10(p), 2) + 3) {
                var round;
                do {
                    var p1 = _subtprob(n, x);
                    var n1 = n + 1;
                    var delta = (p1 - p)
                        / Math.exp((n1 * Math.log(n1 / (n + x * x))
                            + Math.log(n/n1/2/Math.PI) - 1
                            + (1/n1 - 1/n) / 6) / 2);
                    x += delta;
                    round = delta.toFixed(Math.abs(integer(Math.log10(Math.abs(x))-4)));
                } while ((x) && (round != 0));
            }
            return x;
}

function _subu (p) {
    var y = -Math.log(4 * p * (1 - p));
    var x = Math.sqrt(
        y * (1.570796288
        + y * (.03706987906
            + y * (-.8364353589E-3
                + y *(-.2250947176E-3
                    + y * (.6841218299E-5
                        + y * (0.5824238515E-5
                            + y * (-.104527497E-5
                                + y * (.8360937017E-7
                                    + y * (-.3231081277E-8
                                        + y * (.3657763036E-10
                                            + y *.6936233982E-12)))))))))));
    if (p>.5)
        x = -x;
    return x;
}

function _subtprob (n, x) {
    var a;
    var b;
    var w = Math.atan2(x / Math.sqrt(n), 1);
    var z = Math.pow(Math.cos(w), 2);
    var y = 1;
    for (var i = n-2; i >= 2; i -= 2) {
        y = 1 + (i-1) / i * z * y;
    }

    if (n % 2 == 0) {
        a = Math.sin(w)/2;
        b = .5;
    } else {
        a = (n == 1) ? 0 : Math.sin(w)*Math.cos(w)/Math.PI;
        b= .5 + w/Math.PI;
    }
    return max(0, 1 - b - a * y);
}



function stuT(a,b){a=Math.abs(a);var c=a/Math.sqrt(b),d=Math.atan(c);if(1==b)return d/(Math.PI*2);var e=Math.sin(d),f=Math.cos(d);return alpha=b%2==1?1-(d+e*f*stuComp(f*f,2,b-3,-1))/(Math.PI*2):1-e*stuComp(f*f,1,b-3,-1),1-alpha}
function stuComp(a,b,c,d){for(var e=1,f=e,g=b;g<=c;)e=e*a*g/(g-d),f+=e,g+=2;return f}
function pStuT(a,b){for(var c=.5,d=.5,e=0;d>1e-6;)e=1/c-1,d/=2,qt=1-stuT(e,b),qt>a?c-=d:c+=d;return e}
function resConvert(a){var b;return b=a>=0?a+5e-4:a-5e-4}
function easyRoundOf(a,b){if(isNaN(a))return 0;b=Math.pow(10,parseConv(b));var c=Math.round(parseConv(a)*b)/b;return isNaN(c)?0:c}
function parseConv(a){return parseFloat(a)}


function doubleFactorial(n) {
    if(n <= 0) return 1;
    var res = n;
    var tempN = n;
    while (tempN-2 > 0) {
        tempN -= 2;
        res *= tempN;
    }

    return res;
}

function Factorial(n) {
    if(n <= 0) return 1;
    var res = n;
    var tempN = n;
    while (tempN - 1 > 0) {
        tempN -= 1;
        res *= tempN;
    }
    return res;
}

function GammaFunction(N) {
    if (N - N.toFixed(0) == 0) {
        return Factorial(N);
    } else if (Math.abs(N - N.toFixed(0)) == 0.5) {
        return GammaFunctionForHalfes(N);
    }
}

function GammaFunctionForHalfes(N) {
    n = N - 0.5;
    return doubleFactorial(2*n-1)* Math.sqrt(Math.PI)/Math.pow(2, n);
}

function tValueDensity(n, p){
    return GammaFunctionForHalfes((n+1)/2)/(Math.sqrt(n*Math.PI) * GammaFunctionForHalfes(n/2)) * Math.pow(1 + Math.pow(p, 2)/n, -(n+1)/2);
}

function BetaFunction(x, y){
    return GammaFunction(x) * GammaFunction(y) / GammaFunction(x + y);
}

function IncompleteBetaFunction(z, a, b) {
    return Math.pow(z, a)/a * HypergeometricFunction(a, 1 - b, a + 1, z);
}
function RegularizedIncompleteBetaFunction(z, a, b) {
    return IncompleteBetaFunction(z, a, b)/BetaFunction(a, b);
}

function PochhammerSymbol(x, n) {
    return gammaNumber(x + n) / gammaNumber(x);
}

function HypergeometricFunction(a, b, c, z) {
    const iterationsNum = 80;

    var res = 0;

    for (let n = 0; n < iterationsNum; n++) {
        res += (PochhammerSymbol(a, n) * PochhammerSymbol(b, n)/PochhammerSymbol(c, n)) * Math.pow(z, n) / Factorial(n);
    }

    return res;
}

/*function EulerGammaFunction(z) {
    const iterationsNum = 2;
    var res = 1/z;

    for (let n = 0; n < iterationsNum; n++) {
        console.log(res);
        res *= Math.pow(1 + 1/n, z) * Math.pow(1 + z/n, -1);
    }

    return res;
}*/

function product (i, n) { //factorial for integers
    if (n < i) {
      return 1
    }
  
    if (n === i) {
      return n
    }
  
    const half = (n + i) >> 1 // divide (n + i) by 2 and truncate to integer
    return product(i, half) * product(half + 1, n)
}

function gammaNumber (n) {
    let x
  
    if (n - n.toFixed == 0) {
      if (n <= 0) {
        return isFinite(n) ? Infinity : NaN
      }
  
      if (n > 171) {
        return Infinity // Will overflow
      }
  
      return product(1, n - 1)
    }
  
    if (n < 0.5) {
      return Math.PI / (Math.sin(Math.PI * n) * gammaNumber(1 - n))
    }
  
    if (n >= 171.35) {
      return Infinity // will overflow
    }
  
    if (n > 85.0) { // Extended Stirling Approx
      const twoN = n * n
      const threeN = twoN * n
      const fourN = threeN * n
      const fiveN = fourN * n
      return Math.sqrt(2 * Math.PI / n) * Math.pow((n / Math.E), n) *
        (1 + 1 / (12 * n) + 1 / (288 * twoN) - 139 / (51840 * threeN) -
          571 / (2488320 * fourN) + 163879 / (209018880 * fiveN) +
          5246819 / (75246796800 * fiveN * n))
    }
  
    --n
    x = gammaP[0]
    for (let i = 1; i < gammaP.length; ++i) {
      x += gammaP[i] / (n + i)
    }
  
    const t = n + gammaG + 0.5
    return Math.sqrt(2 * Math.PI) * Math.pow(t, n + 0.5) * Math.exp(-t) * x
  }
  
  // TODO: comment on the variables g and p
  
const gammaG = 4.7421875
  
const gammaP = [
    0.99999999999999709182,
    57.156235665862923517,
    -59.597960355475491248,
    14.136097974741747174,
    -0.49191381609762019978,
    0.33994649984811888699e-4,
    0.46523628927048575665e-4,
    -0.98374475304879564677e-4,
    0.15808870322491248884e-3,
    -0.21026444172410488319e-3,
    0.21743961811521264320e-3,
    -0.16431810653676389022e-3,
    0.84418223983852743293e-4,
    -0.26190838401581408670e-4,
    0.36899182659531622704e-5
  ]
  
  // lgamma implementation ref: https://mrob.com/pub/ries/lanczos-gamma.html#code
  
  // log(2 * pi) / 2
const lnSqrt2PI = 0.91893853320467274178
  
const lgammaG = 5 // Lanczos parameter "g"
const lgammaN = 7 // Range of coefficients "n"
  
const lgammaSeries = [
    1.000000000190015,
    76.18009172947146,
    -86.50532032941677,
    24.01409824083091,
    -1.231739572450155,
    0.1208650973866179e-2,
    -0.5395239384953e-5
  ]
  
function lgammaNumber (n) {
    if (n < 0) return NaN
    if (n === 0) return Infinity
    if (!isFinite(n)) return n
  
    if (n < 0.5) {
      // Use Euler's reflection formula:
      // gamma(z) = PI / (sin(PI * z) * gamma(1 - z))
      return Math.log(Math.PI / Math.sin(Math.PI * n)) - lgammaNumber(1 - n)
    }
  
    // Compute the logarithm of the Gamma function using the Lanczos method
  
    n = n - 1
    const base = n + lgammaG + 0.5 // Base of the Lanczos exponential
    let sum = lgammaSeries[0]
  
    // We start with the terms that have the smallest coefficients and largest denominator
    for (let i = lgammaN - 1; i >= 1; i--) {
      sum += lgammaSeries[i] / (n + i)
    }
  
    return lnSqrt2PI + (n + 0.5) * Math.log(base) - base + Math.log(sum)
  }


function tValue(n, t) {
    return 0.5 + t*GammaFunctionForHalfes(0.5*(n+1))/(Math.sqrt(Math.PI * n) * GammaFunctionForHalfes(n/2)) * HypergeometricFunction(0.5, 0.5*(n + 1), 1.5, - Math.pow(t, 2)/n);
}

function tValueI(n, t) {
    return 1 - 0.5*RegularizedIncompleteBetaFunction(n/(Math.pow(t, 2) + n), n/2, 0.5);
}


function StudentsTQuantile(n = 10, p = 0.95) {
    if(n < 1 || p <= 0 || p > 1){
        console.error("Error: invalid n or p.");
        return 0;
    } else if(n == 1){
        return Math.cos(p*Math.PI/2)/Math.sin(p*Math.PI/2);
    }else if(n == 2){
        return Math.sqrt(2 / (p * (2-p)) - 2);
    }

    var a = 1/(n-0.5);
    var b = 48/Math.pow(a, 2);
    var c =((20700*a/b - 98) * a - 16) * a + 96.36;
    var d = ((94.5/(b+c) - 3)/b + 1) * Math.sqrt(a*Math.PI/2) * n;
    var x = d * p;
    var y = Math.pow(x, 2/n);

    if(y > 0.05 + a){
        x = normdev(p * 0.5);
        y = x**2;

        if (n < 5) c = c + 0.3*(n-4.5)*(x+0.6);

        c = (((0.05*d*x - 5) * x - 7) * x - 2) * x + b + c;
        y = (((((0.4*y + 6.3) * y + 36) * y + 94.5)/c - y - 3)/b + 1) * x;
        y = a*(y**2);
        y = (y > 0.002) ? (Math.exp(y) - 1) : (0.5 * (y**2) + y);
    } else{
        y = ((1/(((n + 6)/(n*y) - 0.089*d - 0.822) * (n+2)*3) + 0.5/(n+4)) * y - 1) * (n+1)/(n+2) + 1/y;
    }

    return Math.sqrt(n * y);
}

function normdev(n){
    return Math.exp(-(n**2)/2)/Math.sqrt(2 * Math.PI);
}

function test(n, p){
    console.log(StudentsTQuantile(n, p));
    console.log(tValueDensity(n, p));
    console.log(tValue(n, p));
    console.log(tValueI(n, p));
    console.log(tdistr(n, p));
    console.log(stuT(n, p), stuT(p, n));
}