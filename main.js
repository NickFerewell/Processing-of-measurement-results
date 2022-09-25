function CreateTable(significanceLevel, ...measurements/*samples*/) { //40.15 40.20 40.10 40.15 40.15 40.15 40.10 40.25 40.10 40.15
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
        deviations.push(deviation);
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
    var d2 = Math.sqrt(0.05**2 + d1**2);

    var table = "";
    table += "№" + "  │  " + "x, у.е." + "  │  " + "(x-〈x〉), у.е." + "  │  " + "(x-〈x〉)², у.е." + "\n";
    for (let i = 0; i < measurements.length; i++) {
        table += i + " │ " + measurements[i] + " │ " +  deviations[i] + " │ " +  squaredDeviations[i].toFixed(deviations[i].toString().length - 1 - deviations[i].toString().indexOf(".")) + "\n";
    }
    table += "Σ │ " + measurementsSum + " │ " + deviationsSum + " │ " + squaredDeviationsSum + "\n";
    table += "〈x〉 │ " + meanValue + "\n"; // ̅
    console.log(table)

    answer = "Ответ: (" + meanValue + "±" + d2.toFixed(meanValue.toString().length - 1 - meanValue.toString().indexOf(".")) + ") у.е., n = " + n + ", p = " + significanceLevel + ".";

    console.log(answer);
    table += answer;
    
    return table;
}

function updateTable() { //display, calculate, update
    const significanceLevel = document.getElementById("significanceInput").value;
    const measurements = document.getElementById("measurementsInput").value.split(" ").map(function (element) {
        return +element;
    });
    console.log(significanceLevel, measurements)
    document.getElementById("table").innerText = CreateTable(significanceLevel, ...measurements);
}