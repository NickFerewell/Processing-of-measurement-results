function CreateTableOld(significanceLevel, sistematicError, ...measurements/*samples*/) { //40.15 40.20 40.10 40.15 40.15 40.15 40.10 40.25 40.10 40.15
    //df - degrees of freedom, number of measurements
    n = measurements.length;
    var meanValue = 0;
    measurements.forEach(measurement => {
        meanValue += measurement;
    });
    var measurementsSum = meanValue;
    meanValue /= n;

    var deviations = [];
    var deviationsSum = 0;
    measurements.forEach(measurement => {
        var deviation = measurement - meanValue;
        deviations.push(floorAtIndex(deviation, 4));
        deviationsSum += deviation;
    });

    var squaredDeviations = [];
    var squaredDeviationsSum = 0;
    deviations.forEach(deviation => {
        var squaredDeviation = deviation**2
        squaredDeviations.push(squaredDeviation);
        squaredDeviationsSum += squaredDeviation;
    });

    var Sx = Math.sqrt(squaredDeviationsSum/(n - 1));
    var S = Sx/Math.sqrt(n);
    /*// var t = tdistr(n, significanceLevel)
    var μ = meanValue;
    var x̄ = meanValue*(n/(n-1))
    t = ( x̄ - μ) / (Sx / Math.sqrt(n));*/
    var t;
    if(significanceLevel == 0.95 && n == 10){
        t = 2.262155968610544;
    } else{
        t = tdistr(n-1, (1-significanceLevel)/2);
    }

    var d1 = t * S;
    var d2 = Math.sqrt((sistematicError || 0.05)**2 + d1**2);

    var table = "";
    table += "№" + "  │  " + "x, у.е." + "  │  " + "(x-〈x〉), у.е." + "  │  " + "(x-〈x〉)², у.е." + "\n";
    for (let i = 0; i < measurements.length; i++) {
        table += i+1 + " │ " + measurements[i].toFixed(maxDecimalDigits(measurements)) + " │ " +  deviations[i].toFixed(maxDecimalDigits(deviations)) + " │ " +  squaredDeviations[i].toFixed(maxDecimalDigits(deviations)*2) + "\n";
    }
    table += "Σ │ " + measurementsSum.toFixed(maxDecimalDigits(measurements)) + " │ " + deviationsSum.toFixed(maxDecimalDigits(deviations)) + " │ " + squaredDeviationsSum.toFixed(maxDecimalDigits(deviations)*2) + "\n";
    table += "〈x〉 │ " + meanValue + "\n"; // ̅
    console.log(table)

    table += "S(СКО) - " + S + "\n";
    table += "Погрешность - " + d2 + "\n";

    answer = "Ответ: (" + meanValue + "±" + d2.toFixed(meanValue.toString().length - 1 - meanValue.toString().indexOf(".")) + ") у.е., n = " + n + ", p = " + significanceLevel + ".";

    console.log(answer);
    table += answer;
    
    return table;
}

function CreateTable(significanceLevel, sistematicError, ...measurements/*samples*/) { //40.15 40.20 40.10 40.15 40.15 40.15 40.10 40.25 40.10 40.15
    //df - degrees of freedom, number of measurements
    n = measurements.length;
    var meanValue = 0;
    measurements.forEach(measurement => {
        meanValue += measurement;
    });
    var measurementsSum = meanValue;
    meanValue /= n;

    var deviations = [];
    var deviationsSum = 0;
    measurements.forEach(measurement => {
        var deviation = measurement - meanValue;
        deviations.push(floorAtIndex(deviation, 4));
        deviationsSum += deviation;
    });

    var squaredDeviations = [];
    var squaredDeviationsSum = 0;
    deviations.forEach(deviation => {
        var squaredDeviation = deviation**2
        squaredDeviations.push(squaredDeviation);
        squaredDeviationsSum += squaredDeviation;
    });

    var Sx = Math.sqrt(squaredDeviationsSum/(n - 1));
    var S = Sx/Math.sqrt(n);
    var t;
    if(significanceLevel == 0.95 && n == 10){
        t = 2.262155968610544;
    } else{
        t = tdistr(n-1, (1-significanceLevel)/2);
    }

    var d1 = t * S;
    var d2 = Math.sqrt(sistematicError**2 + d1**2);

    var Mtable = [];
    Mtable.push(["№", "x, у.е.","(x-〈x〉), у.е.","(x-〈x〉)², у.е."]);
    console.log(Mtable)
    for (let i = 0; i < measurements.length; i++) {
        Mtable.push([i+1, measurements[i].toFixed(maxDecimalDigits(measurements)), deviations[i].toFixed(maxDecimalDigits(deviations)), squaredDeviations[i].toFixed(maxDecimalDigits(deviations)*2)]);
    }
    console.log(Mtable)
    Mtable.push(["Σ", measurementsSum.toFixed(maxDecimalDigits(measurements)), deviationsSum.toFixed(maxDecimalDigits(deviations)), squaredDeviationsSum.toFixed(maxDecimalDigits(deviations)*2)]);
    Mtable.push(["〈x〉", meanValue]); // ̅
    console.log(Mtable)
    var tableRes = table(Mtable, {hsep: "│"}) + "\n";

    tableRes += "S(СКО) - " + S + "\n";
    tableRes += "Погрешность - " + d2 + "\n";

    answer = "Ответ: (" + meanValue + "±" + d2.toFixed(meanValue.toString().length - 1 - meanValue.toString().indexOf(".")) + ") у.е., n = " + n + ", p = " + significanceLevel + ".";

    console.log(answer);
    tableRes += answer;
    
    return tableRes;
}

function updateTable() { //display, calculate, update
    const significanceLevel = document.getElementById("significanceInput").value || +document.getElementById("significanceInput").getAttribute("placeholder");
    const sistematicError = document.getElementById("sistematicErrorInput").value || +document.getElementById("sistematicErrorInput").getAttribute("placeholder");
    const measurements = document.getElementById("measurementsInput").value.split(" ").map(function (element) {
        return +element;
    });
    console.log(significanceLevel, sistematicError, measurements)
    document.getElementById("table").innerText = CreateTable(significanceLevel, sistematicError, ...measurements);
}

function floorAtIndex(n, i) {
    return Math.floor(n * (10**i))/(10**i);
}

function maxDecimalDigits(arr){
    max = 0;
    arr.forEach((n) => {
        max = Math.max(max, n.toString().length - 1 - n.toString().indexOf("."));
    })
    return max;
}