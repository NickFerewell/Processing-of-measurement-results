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
                    round = delta.toFixed(Math.abs((Math.log10(Math.abs(x))-4).toFixed(0)));
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
    return Math.max(0, 1 - b - a * y);
}