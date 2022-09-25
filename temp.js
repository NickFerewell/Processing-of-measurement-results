
        var SIGNIFICANT = 5; // number of significant digits to be returned

        function chisqrdistr ($n, $p) {
            if ($n <= 0 || Math.abs($n) - Math.abs(integer($n)) != 0) {
                throw("Invalid n: $n\n"); /* degree of freedom */
            }
            if ($p <= 0 || $p > 1) {
                throw("Invalid p: $p\n");
            }
            return precision_string(_subchisqr($n-0, $p-0));
        }

        function udistr ($p) {
            if ($p > 1 || $p <= 0) {
                throw("Invalid p: $p\n");
            }
            return precision_string(_subu($p-0));
        }

        function tdistr ($n, $p) {
            if ($n <= 0 || Math.abs($n) - Math.abs(integer($n)) != 0) {
                throw("Invalid n: $n\n");
            }
            if ($p <= 0 || $p >= 1) {
                throw("Invalid p: $p\n");
            }
            return precision_string(_subt($n-0, $p-0));
        }

        function fdistr ($n, $m, $p) {
            if (($n<=0) || ((Math.abs($n)-(Math.abs(integer($n))))!=0)) {
                throw("Invalid n: $n\n"); /* first degree of freedom */
            }
            if (($m<=0) || ((Math.abs($m)-(Math.abs(integer($m))))!=0)) {
                throw("Invalid m: $m\n"); /* second degree of freedom */
            }
            if (($p<=0) || ($p>1)) {
                throw("Invalid p: $p\n");
            }
            return precision_string(_subf($n-0, $m-0, $p-0));
        }

        function uprob ($x) {
            return precision_string(_subuprob($x-0));
        }

        function chisqrprob ($n,$x) {
            if (($n <= 0) || ((Math.abs($n) - (Math.abs(integer($n)))) != 0)) {
                throw("Invalid n: $n\n"); /* degree of freedom */
            }
            return precision_string(_subchisqrprob($n-0, $x-0));
        }

        function tprob ($n, $x) {
            if (($n <= 0) || ((Math.abs($n) - Math.abs(integer($n))) !=0)) {
                throw("Invalid n: $n\n"); /* degree of freedom */
            }
            return precision_string(_subtprob($n-0, $x-0));
        }

        function fprob ($n, $m, $x) {
            if (($n<=0) || ((Math.abs($n)-(Math.abs(integer($n))))!=0)) {
                throw("Invalid n: $n\n"); /* first degree of freedom */
            }
            if (($m<=0) || ((Math.abs($m)-(Math.abs(integer($m))))!=0)) {
                throw("Invalid m: $m\n"); /* second degree of freedom */
            }
            return precision_string(_subfprob($n-0, $m-0, $x-0));
        }


        function _subfprob ($n, $m, $x) {
            var $p;

            if ($x<=0) {
                $p=1;
            } else if ($m % 2 == 0) {
                var $z = $m / ($m + $n * $x);
                var $a = 1;
                for (var $i = $m - 2; $i >= 2; $i -= 2) {
                    $a = 1 + ($n + $i - 2) / $i * $z * $a;
                }
                $p = 1 - Math.pow((1 - $z), ($n / 2) * $a);
            } else if ($n % 2 == 0) {
                var $z = $n * $x / ($m + $n * $x);
                var $a = 1;
                for (var $i = $n - 2; $i >= 2; $i -= 2) {
                    $a = 1 + ($m + $i - 2) / $i * $z * $a;
                }
                $p = Math.pow((1 - $z), ($m / 2)) * $a;
            } else {
                var $y = Math.atan2(Math.sqrt($n * $x / $m), 1);
                var $z = Math.pow(Math.sin($y), 2);
                var $a = ($n == 1) ? 0 : 1;
                for (var $i = $n - 2; $i >= 3; $i -= 2) {
                    $a = 1 + ($m + $i - 2) / $i * $z * $a;
                }
                var $b = Math.PI;
                for (var $i = 2; $i <= $m - 1; $i += 2) {
                    $b *= ($i - 1) / $i;
                }
                var $p1 = 2 / $b * Math.sin($y) * Math.pow(Math.cos($y), $m) * $a;

                $z = Math.pow(Math.cos($y), 2);
                $a = ($m == 1) ? 0 : 1;
                for (var $i = $m-2; $i >= 3; $i -= 2) {
                    $a = 1 + ($i - 1) / $i * $z * $a;
                }
                $p = max(0, $p1 + 1 - 2 * $y / Math.PI
                    - 2 / Math.PI * Math.sin($y) * Math.cos($y) * $a);
            }
            return $p;
        }


        function _subchisqrprob ($n,$x) {
            var $p;

            if ($x <= 0) {
                $p = 1;
            } else if ($n > 100) {
                $p = _subuprob((Math.pow(($x / $n), 1/3)
                    - (1 - 2/9/$n)) / Math.sqrt(2/9/$n));
            } else if ($x > 400) {
                $p = 0;
            } else {
                var $a;
                var $i;
                var $i1;
                if (($n % 2) != 0) {
                    $p = 2 * _subuprob(Math.sqrt($x));
                    $a = Math.sqrt(2/Math.PI) * Math.exp(-$x/2) / Math.sqrt($x);
                    $i1 = 1;
                } else {
                    $p = $a = Math.exp(-$x/2);
                    $i1 = 2;
                }

                for ($i = $i1; $i <= ($n-2); $i += 2) {
                    $a *= $x / $i;
                    $p += $a;
                }
            }
            return $p;
        }

        function _subu ($p) {
            var $y = -Math.log(4 * $p * (1 - $p));
            var $x = Math.sqrt(
                $y * (1.570796288
                + $y * (.03706987906
                    + $y * (-.8364353589E-3
                        + $y *(-.2250947176E-3
                            + $y * (.6841218299E-5
                                + $y * (0.5824238515E-5
                                    + $y * (-.104527497E-5
                                        + $y * (.8360937017E-7
                                            + $y * (-.3231081277E-8
                                                + $y * (.3657763036E-10
                                                    + $y *.6936233982E-12)))))))))));
            if ($p>.5)
                $x = -$x;
            return $x;
        }

        function _subuprob ($x) {
            var $p = 0; /* if ($absx > 100) */
            var $absx = Math.abs($x);

            if ($absx < 1.9) {
                $p = Math.pow((1 +
                    $absx * (.049867347
                        + $absx * (.0211410061
                            + $absx * (.0032776263
                                + $absx * (.0000380036
                                    + $absx * (.0000488906
                                        + $absx * .000005383)))))), -16)/2;
            } else if ($absx <= 100) {
                for (var $i = 18; $i >= 1; $i--) {
                    $p = $i / ($absx + $p);
                }
                $p = Math.exp(-.5 * $absx * $absx)
                    / Math.sqrt(2 * Math.PI) / ($absx + $p);
            }

            if ($x<0)
                $p = 1 - $p;
            return $p;
        }


        function _subt ($n, $p) {

            if ($p >= 1 || $p <= 0) {
                throw("Invalid p: $p\n");
            }

            if ($p == 0.5) {
                return 0;
            } else if ($p < 0.5) {
                return - _subt($n, 1 - $p);
            }

            var $u = _subu($p);
            var $u2 = Math.pow($u, 2);

            var $a = ($u2 + 1) / 4;
            var $b = ((5 * $u2 + 16) * $u2 + 3) / 96;
            var $c = (((3 * $u2 + 19) * $u2 + 17) * $u2 - 15) / 384;
            var $d = ((((79 * $u2 + 776) * $u2 + 1482) * $u2 - 1920) * $u2 - 945)
                / 92160;
            var $e = (((((27 * $u2 + 339) * $u2 + 930) * $u2 - 1782) * $u2 - 765) * $u2
                + 17955) / 368640;

            var $x = $u * (1 + ($a + ($b + ($c + ($d + $e / $n) / $n) / $n) / $n) / $n);

            if ($n <= Math.pow(log10($p), 2) + 3) {
                var $round;
                do {
                    var $p1 = _subtprob($n, $x);
                    var $n1 = $n + 1;
                    var $delta = ($p1 - $p)
                        / Math.exp(($n1 * Math.log($n1 / ($n + $x * $x))
                            + Math.log($n/$n1/2/Math.PI) - 1
                            + (1/$n1 - 1/$n) / 6) / 2);
                    $x += $delta;
                    $round = round_to_precision($delta, Math.abs(integer(log10(Math.abs($x))-4)));
                } while (($x) && ($round != 0));
            }
            return $x;
        }
        function _subtprob ($n, $x) {
            var $a;
            var $b;
            var $w = Math.atan2($x / Math.sqrt($n), 1);
            var $z = Math.pow(Math.cos($w), 2);
            var $y = 1;
            for (var $i = $n-2; $i >= 2; $i -= 2) {
                $y = 1 + ($i-1) / $i * $z * $y;
            }

            if ($n % 2 == 0) {
                $a = Math.sin($w)/2;
                $b = .5;
            } else {
                $a = ($n == 1) ? 0 : Math.sin($w)*Math.cos($w)/Math.PI;
                $b= .5 + $w/Math.PI;
            }
            return max(0, 1 - $b - $a * $y);
        }

        function _subf ($n, $m, $p) {
            var $x;

            if ($p >= 1 || $p <= 0) {
                throw("Invalid p: $p\n");
            }

            if ($p == 1) {
                $x = 0;
            } else if ($m == 1) {
                $x = 1 / Math.pow(_subt($n, 0.5 - $p / 2), 2);
            } else if ($n == 1) {
                $x = Math.pow(_subt($m, $p/2), 2);
            } else if ($m == 2) {
                var $u = _subchisqr($m, 1 - $p);
                var $a = $m - 2;
                $x = 1 / ($u / $m * (1 +
                    (($u - $a) / 2 +
                        (((4 * $u - 11 * $a) * $u + $a * (7 * $m - 10)) / 24 +
                            (((2 * $u - 10 * $a) * $u + $a * (17 * $m - 26)) * $u
                                - $a * $a * (9 * $m - 6)
                            )/48/$n
                        )/$n
                    )/$n));
            } else if ($n > $m) {
                $x = 1 / _subf2($m, $n, 1 - $p)
            } else {
                $x = _subf2($n, $m, $p)
            }
            return $x;
        …