function test() {
	x = document.getElementById("totalNumMonomers").value
	alert(x);
}
function simulate() {
	//Get all relevant inputs
	var numUniqueMonomers = 2
	var totalNumMonomers = parseInt(document.getElementById("totalNumMonomers").value);
	var mRatio = parseInt(document.getElementById("mRatio").value);
	var conversion = parseInt(document.getElementById("conversion").value);
	var numRowsToShow = parseInt(document.getElementById("numRowsToShow").value);
	var graph1TypeObj = document.getElementById("graph1Type");
	var graph2TypeObj = document.getElementById("graph2Type");
	var graph1Type = graph1TypeObj.options[graph1TypeObj.selectedIndex].value;
	var graph2Type = graph2TypeObj.options[graph2TypeObj.selectedIndex].value;
	var hist1Monomer = parseInt(document.getElementById("hist1Monomer").value);
	var hist2Monomer = parseInt(document.getElementById("hist2Monomer").value);
	var monomer1Ratio = parseFloat(document.getElementById("monomer1Ratio").value);
	var monomer2Ratio = parseFloat(document.getElementById("monomer2Ratio").value);
	var monomer1RR = parseFloat(document.getElementById("monomer1RR").value);
	var monomer2RR = parseFloat(document.getElementById("monomer2RR").value);
	//Initial variable calculations
	console.log("Type of monomer1Ratio: ", typeof(monomer1Ratio));
	var polymerLength = Math.floor(mRatio * conversion / 100);
	var numPolymers = Math.floor(totalNumMonomers / mRatio);
	var monomerRatioList = [monomer1Ratio, monomer2Ratio];
	var rrList = [[monomer1RR, 1], [1, monomer2RR]];
	var monomerAmountsList = getMonomerAmounts(monomerRatioList, totalNumMonomers);
	console.log(monomerAmountsList);
	//Initiate chains, all polymer chains are represented as arrays and stored in another array, polymerArray
	var polymerArray = [];
	var currNumPolymers;
	for (currNumPolymers = 0; currNumPolymers < numPolymers; currNumPolymers++) {
		var initChoices = [];
		var initWeightList = [];
		//if it is a 2-monomer system, use the Mayo Lewis Equation to Determine Starting Monomers
		//Weights that determine probabilty of monomer starting the chain are calculated
		if (numUniqueMonomers == 2) {
			var f1 = monomerAmountsList[0];
			var f2 = monomerAmountsList[1];
			var r1 = rrList[0][0];
			var r2 = rrList[1][1];
			var initWeight = (r1*f1**2 + f1*f2) / (r1*f1**2 + 2*f1*f2 + r2*f2**2);
			initChoices.push(1);
			initChoices.push(2);
			initWeightList.push(initWeight);
			initWeightList.push(1 - initWeight);
		} 
		//Use a weighted random selector to choose the starting monomer
		var startingMonomer = parseInt(weightedRand(initChoices, initWeightList));
		//console.log("typeof startingMonomer: ", typeof(startingMonomer));
		//console.log("startingMonomer: ", startingMonomer);
		//Create a new polymer chain with the selected starting monomer
		polymerArray.push([startingMonomer]);
		//Remove one of that monomer from the pool
		monomerAmountsList[startingMonomer - 1] --;
	}
	//console.log("polymerArray: ", polymerArray);
	//Propogate chains
	//First loop: propagate until the polymer length reaches expected value
	for (var currLength = 1; currLength < polymerLength; currLength++) {
		//Second loop: during each cycle of propagation, add one monomer to each of the polymer chains
		for (var i = 0; i < polymerArray.length; i++) {
			//console.log("polymerArray: ", polymerArray);
			var polymer = polymerArray[i];
			var propagationChoices = [];
			var propagationWeightList = [];
			//Third loop: cycle through each unique monomer and calculate probability weight of the monomer to be added to 
			// the growing polymer chain
			for (var monomerID = 1; monomerID < numUniqueMonomers + 1; monomerID++) {
				//console.log("rrList: ", rrList);
				//console.log("polymer: ", polymer);
				//console.log("last monomer: ", polymer[polymer.length - 1] - 1);
				var reactivityRatio = rrList[polymer[polymer.length - 1] - 1][monomerID - 1];
				var propWeight = monomerAmountsList[monomerID - 1] * reactivityRatio;
				//console.log("monomerID: ",monomerID, "propWeight: ", propWeight);
				propagationChoices.push(monomerID);
				propagationWeightList.push(propWeight);
			}
			//Using weights, randomly select propagating monomer and add it to the polymer chain array
			//console.log("propagationSpec:", propagationSpec);
			var nextMonomer = parseInt(weightedRand(propagationChoices, propagationWeightList));
			polymer.push(nextMonomer);
			//Remove that monomer from the reactant pool
			monomerAmountsList[nextMonomer - 1]--;
		}
	}
	//console.log("final polymerArray: ", polymerArray);
	var chartData = getMonomerComposition(polymerArray, numUniqueMonomers, polymerLength);
	updateChart(chartData);
	console.log("chartData: ", chartData);

}
function updateChart(chartData) {
	chart.dataProvider = chartData;
	chart.validateData();
}
function getMonomerAmounts(monomerRatioList, totalNumMonomers) {
	var monomerAmountsList = [];
	var totalWeight = sumArray(monomerRatioList);
	for (var i = 0; i < monomerRatioList.length; i++) {
		var monomerRatio = monomerRatioList[i];
		var weight = monomerRatio / totalWeight;
		var monomerAmount = Math.ceil(totalNumMonomers * weight);
		//console.log("monomerRatio: ", monomerRatio);
		monomerAmountsList.push(monomerAmount);
	}
	return monomerAmountsList;
}
//Returns chartData formatted information for percent composition of each unique monomer at each index.
function getMonomerComposition(polymerArray, numUniqueMonomers, polymerLength) {
	var fullCompositionList = [];
	for (var monomerID = 1; monomerID <= numUniqueMonomers; monomerID++) {
		var singleCompositionList = [];
		for (var positionIndex = 0; positionIndex < polymerLength; positionIndex++) {
			var monomerCount = 0;
			for (var i = 0; i < polymerArray.length; i++) {
				polymer = polymerArray[i];
				if (polymer[positionIndex] == monomerID) {
					monomerCount++;
				}
			}
			var monomerRatio = monomerCount / polymerArray.length;
			singleCompositionList.push(monomerRatio);
		}
		var xyDataPair = [createRangeArray(polymerLength), singleCompositionList];
		fullCompositionList.push(xyDataPair);
	}
	var chartData = [];
	for (var i = 0; i <= polymerLength; i++) {
		var obj = {};
		for (var monomerID = 1; monomerID <= numUniqueMonomers; monomerID++) {
			var xID = "x" + monomerID;
			var yID = "y" + monomerID;
			console.log("xID: ", xID);
			obj[xID] = fullCompositionList[monomerID - 1][0][i];
			obj[yID] = fullCompositionList[monomerID - 1][1][i];
		}
		chartData.push(obj);
	}
	return chartData;
}
function convertToAmChartData(x, y, id) {
	var chartData = [];
	for( var i = 0; i < x.length; i++ ) {
	  chartData.push( {
	    "x" : x[ i ],
	    "y" : y[ i ]
	  } )
	}
	return chartData;
}
function show() {
  var chartData = [{
    "x1": 1,
    "y1": 0.5,
    "x2": 1,
    "y2": 2.2
  }, {
    "x1": 2,
    "y1": 1.3,
    "x2": 2,
    "y2": 4.9
  }, {
    "x1": 3,
    "y1": 2.3,
    "x2": 3,
    "y2": 5.1
  }, {
    "x1": 4,
    "y1": 2.8,
    "x2": 4,
    "y2": 5.3
  }, {
    "x1": 5,
    "y1": 3.5,
    "x2": 5,
    "y2": 6.1
  }, {
    "x1": 6,
    "y1": 5.1,
    "x2": 6,
    "y2": 8.3
  }, {
    "x1": 7,
    "y1": 6.7,
    "x2": 7,
    "y2": 10.5
  }, {
    "x1": 8,
    "y1": 8,
    "x2": 8,
    "y2": 12.3
  }, {
    "x1": 9,
    "y1": 8.9,
    "x2": 9,
    "y2": 14.5
  }, {
    "x1": 10,
    "y1": 9.7,
    "x2": 10,
    "y2": 15
  }, {
    "x1": 11,
    "y1": 10.4,
    "x2": 11,
    "y2": 18.8
  }, {
    "x1": 12,
    "y1": 11.7,
    "x2": 12,
    "y2": 19
  }];
  console.log("Original chartData: ", chartData);

  // XY CHART
  chart = new AmCharts.AmXYChart();

  chart.dataProvider = chartData;
  chart.startDuration = 1;

  // AXES
  // X
  var xAxis = new AmCharts.ValueAxis();
  xAxis.title = "X Axis";
  xAxis.position = "bottom";
  xAxis.dashLength = 1;
  xAxis.axisAlpha = 0;
  xAxis.autoGridCount = true;
  chart.addValueAxis(xAxis);

  // Y
  var yAxis = new AmCharts.ValueAxis();
  yAxis.position = "left";
  yAxis.title = "Y Axis";
  yAxis.dashLength = 1;
  yAxis.axisAlpha = 0;
  yAxis.autoGridCount = true;
  chart.addValueAxis(yAxis);

  // GRAPHS
  // triangles up
  var graph1 = new AmCharts.AmGraph();
  graph1.lineColor = "#FF6600";
  graph1.balloonText = "x:[[x]] y:[[y]]";
  graph1.xField = "x1";
  graph1.yField = "y1";
  graph1.lineAlpha = 0;
  graph1.bullet = "triangleUp";
  chart.addGraph(graph1);

  // triangles down
  var graph2 = new AmCharts.AmGraph();
  graph2.lineColor = "#FCD202";
  graph2.balloonText = "x:[[x]] y:[[y]]";
  graph2.xField = "x2";
  graph2.yField = "y2";
  graph2.lineAlpha = 0;
  graph2.bullet = "triangleDown";
  chart.addGraph(graph2);

  // first trend line
  var trendLine = new AmCharts.TrendLine();
  trendLine.lineColor = "#FF6600";
  trendLine.initialXValue = 1;
  trendLine.initialValue = 2;
  trendLine.finalXValue = 12;
  trendLine.finalValue = 11;
  chart.addTrendLine(trendLine);

  // second trend line
  trendLine = new AmCharts.TrendLine();
  trendLine.lineColor = "#FCD202";
  trendLine.initialXValue = 1;
  trendLine.initialValue = 1;
  trendLine.finalXValue = 12;
  trendLine.finalValue = 19;
  chart.addTrendLine(trendLine);

  // CURSOR
  var chartCursor = new AmCharts.ChartCursor();
  chart.addChartCursor(chartCursor);

  // SCROLLBAR

  var chartScrollbar = new AmCharts.ChartScrollbar();
  chartScrollbar.scrollbarHeight = 5;
  chartScrollbar.offset = 15
  chart.addChartScrollbar(chartScrollbar);
  // WRITE
  chart.write("chartdiv");
}
//Creates an array from [1, 2, 3,......, N]
function createRangeArray(length) {
	var array = [];
	for (var i = 1; i <= length; i++) {
		array.push(i);
	}
	return array;
}
//Sums the numerical contents of an array
function sumArray(array) {
	return array.reduce(sum, 0)
}
function sum(a, b) {
	return a + b;
}
function rand(min, max) {
    return Math.random() * (max - min) + min;
};
 
function weightedRand(list, weight) {
    var total_weight = weight.reduce(function (prev, cur, i, arr) {
        return prev + cur;
    });
     
    var random_num = rand(0, total_weight);
    var weight_sum = 0;
    //console.log(random_num)
     
    for (var i = 0; i < list.length; i++) {
        weight_sum += weight[i];
        weight_sum = +weight_sum.toFixed(2);
         
        if (random_num <= weight_sum) {
            return list[i];
        }
    }
     
    // end of function
};