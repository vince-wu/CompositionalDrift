function simulate() {
	//Get all relevant inputs
	numUniqueMonomers = 2
	var totalNumMonomers = parseInt(document.getElementById("totalNumMonomers").value);
	var mRatio = parseInt(document.getElementById("mRatio").value);
	var conversion = parseInt(document.getElementById("conversion").value);
	var numRowsToShow = parseInt(document.getElementById("numRowsToShow").value);
	var graphTypeObj = document.getElementById("graph1Type");
	var graphType = graphTypeObj.options[graphTypeObj.selectedIndex].value;
	var monomer1Ratio = parseFloat(document.getElementById("monomer1Ratio").value);
	var monomer2Ratio = parseFloat(document.getElementById("monomer2Ratio").value);
	var monomer1RR = parseFloat(document.getElementById("monomer1RR").value);
	var monomer2RR = parseFloat(document.getElementById("monomer2RR").value);
	//Initial variable calculations
	console.log("Type of monomer1Ratio: ", typeof(monomer1Ratio));
	polymerLength = Math.floor(mRatio * conversion / 100);
	var numPolymers = Math.floor(totalNumMonomers / mRatio);
	var monomerRatioList = [monomer1Ratio, monomer2Ratio];
	var rrList = [[monomer1RR, 1], [1, monomer2RR]];
	var monomerAmountsList = getMonomerAmounts(monomerRatioList, totalNumMonomers);
	initialMonomerAmountList = getMonomerAmounts(monomerRatioList, totalNumMonomers);
	console.log("initialMonomerAmountList1: ", initialMonomerAmountList);
	//console.log(monomerAmountsList);
	//Initiate chains, all polymer chains are represented as arrays and stored in another array, polymerArray
	polymerArray = [];
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
	setGraph(graphType);
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
	chartData = convertToChartData(fullCompositionList, numUniqueMonomers, polymerLength);
	return chartData;
}
function getPercentageMonomer(polymerArray, numUniqueMonomers, polymerLength, initialMonomerAmountList) {
	var fullPercentageList = [];
	for (var monomerID = 1; monomerID <= numUniqueMonomers; monomerID++) {
		var singlePercentageList = [];
		var monomerRemaining = initialMonomerAmountList[monomerID -1];
		for (var positionIndex = 0; positionIndex < polymerLength; positionIndex++) {
			for (var i = 0; i < polymerArray.length; i++) {
				polymer = polymerArray[i];
				if (polymer[positionIndex] == monomerID) {
					monomerRemaining--;
				}
			}
			//console.log("initialMonomerAmountList: ", initialMonomerAmountList);
			var monomerPercentage = monomerRemaining / initialMonomerAmountList[monomerID -1];
			singlePercentageList.push(monomerPercentage);
		}
		var xyDataPair = [createRangeArray(polymerLength), singlePercentageList];
		fullPercentageList.push(xyDataPair);
	}
	chartData = convertToChartData(fullPercentageList, numUniqueMonomers, polymerLength);
	return chartData;
}
function getMonomerSeparation(polymerArray, numUniqueMonomers, polymerLength, monomerID){
	//console.log("monomerID: ", monomerID);
	//console.log("typeof monomerID: ", typeof(monomerID));
	var fullSeparationData = new Array(polymerLength + 1).fill(0)
	var largestBlock = 0;
	for (var i = 0; i < polymerArray.length; i++) {
		polymer = polymerArray[i];
		var numConsecutive = 0;
		for (var positionIndex = 0; positionIndex < polymerLength; positionIndex++) {
			//console.log("numConsecutive: ", numConsecutive);
			currMonomerID = polymer[positionIndex];
			//console.log("currMonomerID: ", currMonomerID);
			if (currMonomerID != monomerID && numConsecutive > 0) {
				fullSeparationData[numConsecutive] += numConsecutive;
				if (numConsecutive > largestBlock) {
					largestBlock = numConsecutive;
				}
				numConsecutive = 0;
			}
			if (currMonomerID == monomerID) {
				numConsecutive += 1;
			}
			if (positionIndex == polymerLength - 1) {
				if (numConsecutive != 0) {
					fullSeparationData[numConsecutive] += numConsecutive;
				}
			}
		}
	}
	console.log("largestBlock: ", largestBlock);
	fullSeparationData = fullSeparationData.slice(0, largestBlock + 2);
	var chartData = [];
	for (var blockSize = 0; blockSize < fullSeparationData.length; blockSize++) {
		var obj = {};
		var xID = "blockSize" + monomerID;
		var xCounts = "counts" + monomerID;
		obj[xID] = blockSize;
		obj[xCounts] = fullSeparationData[blockSize];
		chartData.push(obj);
	}
	return chartData;		
}
function getPolymerCompostion(polymerArray, numUniqueMonomers, polymerLength) {
	var fullCompositionList = [];
	for (var monomerID = 1; monomerID <= numUniqueMonomers; monomerID++) {
		var singleCompositionList = [];
		var totalMonomerCount = 0;
		for (var positionIndex = 0; positionIndex < polymerLength; positionIndex++) {
			for (var i = 0; i < polymerArray.length; i++) {
				polymer = polymerArray[i];
				if (polymer[positionIndex] == monomerID) {
					totalMonomerCount++;
				}
			}
			var polymerComposition = totalMonomerCount / (polymerArray.length * (positionIndex + 1));
			singleCompositionList.push(polymerComposition);
		}
		var xyDataPair = [createRangeArray(polymerLength), singleCompositionList];
		fullCompositionList.push(xyDataPair);
	}
	chartData = convertToChartData(fullCompositionList, numUniqueMonomers, polymerLength);
	return chartData;
}
function convertToChartData(fullList, numUniqueMonomers, polymerLength) {
	chartData = [];
	for (var i = 0; i <= polymerLength; i++) {
		var obj = {};
		for (var monomerID = 1; monomerID <= numUniqueMonomers; monomerID++) {
			var xID = "x" + monomerID;
			var yID = "y" + monomerID;
			//console.log("xID: ", xID);
			obj[xID] = fullList[monomerID - 1][0][i];
			obj[yID] = fullList[monomerID - 1][1][i];
		}
		chartData.push(obj);
	}
	return chartData;
}

function switchCharts(chart) {
	document.getElementById("chartdiv").style.display = "none";
	document.getElementById("chartdiv2").style.display = "block";
	chart.invalidateSize();
	chart.animateAgain();
}
function createHistChart(chartData, monomerID) {
	chart2 = new AmCharts.AmSerialChart();
	var obj = chartData[0];
	keyList = Object.keys(obj);
	console.log("keyList: ", keyList);
	//Histogram
	chart2.dataProvider = chartData;
	chart2.categoryField = keyList[0];
	chart2.startDuration = 1;
	chart2.sequencedAnimation = false; 

	//Value Axis
	var valueAxis = new AmCharts.ValueAxis();
	valueAxis.position = "bottom";
	valueAxis.minimum = 0;
	valueAxis.title = "Normalized Separation";
	chart2.addValueAxis(valueAxis);

	//X Axis

	chart2.categoryAxis.title = "Monomer " + monomerID + " Block Size";

	//LEGEND

	var legend = new AmCharts.AmLegend();
	legend.position = "bottom";
	legend.align = "center";
	legend.markerType = "square";
	legend.useGraphSettings = true;
	//chart2.addLegend(legend);

	//Graph
	var graph = new AmCharts.AmGraph();
	graph.type = "column";
	graph.fillAlphas = 1;
	graph.valueField = keyList[1];
	graph.balloonText = "[[value]]";
	chart2.addGraph(graph);
	chart2.write("chartdiv2");
}
function createXYGraph() {
	colorArray = ["#FF6600", "#FCD202", "#B0DE09", "#0D8ECF", "#2A0CD0", 
	"#CD0D74", "#CC0000", "#00CC00", "#0000CC", "#DDDDDD", "#999999", "#333333", "#990000"];
	var chartData = getMonomerComposition(polymerArray, numUniqueMonomers, polymerLength);
	// XY CHART
	chart = new AmCharts.AmXYChart();
	chart.dataProvider = chartData;
	chart.startDuration = 0;
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

	//Remove all old graphs
	var graphList = chart.graphs;
	for (var i = 0; i < graphList.length; i++) {
		chart.removeGraph(graphList[i]);
	}
	//Add graphs depending on number of unique monomers
	for (var graphNumber = 1; graphNumber <= numUniqueMonomers; graphNumber++) {
		var graph = new AmCharts.AmGraph();
	  	graph.lineColor = colorArray[graphNumber - 1];
		graph.balloonText = "x:[[x]] y:[[y]]";
		graph.xField = "x" + graphNumber;
		graph.yField = "y" + graphNumber;
		graph.lineAlpha = 0;
		graph.bullet = "triangleDown";
		chart.addGraph(graph);
	}
	//LEGEND
	var legend = new AmCharts.AmLegend();
	legend.position = "right";
	legend.align = "center";
	legend.markerType = "square";
	legend.useGraphSettings = true;
	chart.addLegend(legend);

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
function setGraph(type) {
	//console.log("type: ", type);
	createXYGraph();
	var chartData;
	switch (type) {
		case "Monomer Occurences":
			chartData = getMonomerComposition(polymerArray, numUniqueMonomers, polymerLength);
			chart.dataProvider = chartData;
			chart.valueAxes[0].title = "Monomer Index";
			chart.valueAxes[1].title = "Monomer Occurence";
			chart.valueAxes[0].minimum = 0;
			chart.valueAxes[0].maximum = polymerLength;
			chart.valueAxes[1].minimum = 0;
			chart.valueAxes[1].maximum = 1;
			for (var i = 0; i < numUniqueMonomers; i++) {
				chart.graphs[i].title = "Monomer " + (i + 1);
			}
			chart.validateData();
			break;
		case "Percentage Monomer":
			chartData = getPercentageMonomer(polymerArray, numUniqueMonomers, polymerLength, initialMonomerAmountList);
			chart.dataProvider = chartData;
			chart.valueAxes[0].title = "Monomer Index";
			chart.valueAxes[1].title = "Unreacted Monomer Remaining";
			chart.valueAxes[0].minimum = 0;
			chart.valueAxes[0].maximum = polymerLength;
			chart.valueAxes[1].minimum = 0;
			chart.valueAxes[1].maximum = 1;
			for (var i = 0; i < numUniqueMonomers; i++) {
				chart.graphs[i].title = "Monomer " + (i + 1);
			}
			chart.validateData();
			break;
		case "Polymer Compositions":
			chartData = getPolymerCompostion(polymerArray, numUniqueMonomers, polymerLength);
			chart.dataProvider = chartData;
			chart.valueAxes[0].title = "Monomer index";
			chart.valueAxes[1].title = "Instantaneous Composition of Polymer Chain";
			chart.valueAxes[0].minimum = 0;
			chart.valueAxes[0].maximum = polymerLength;
			chart.valueAxes[1].minimum = 0;
			chart.valueAxes[1].maximum = 1;
			for (var i = 0; i < numUniqueMonomers; i++) {
				chart.graphs[i].title = "Monomer " + (i + 1);
			}
			chart.validateData();
			break;
		case "Monomer Separation":
			var histMonomer = parseInt(document.getElementById("hist1Monomer").value);
			chartData = getMonomerSeparation(polymerArray, numUniqueMonomers, polymerLength, histMonomer);
			createHistChart(chartData, histMonomer);
			document.getElementById("chartdiv").style.display = "none";
			document.getElementById("chartdiv2").style.display = "block";
			chart2.invalidateSize();
			chart2.animateAgain();
			return;
	}
	document.getElementById("chartdiv2").style.display = "none";
	document.getElementById("chartdiv").style.display = "block";
	console.log("got here");
	chart.invalidateSize();
}
function setHistGraph(monomerID) {
	monomerID = parseInt(monomerID);
	console.log("monomerID: ", monomerID);
	var chartData = getMonomerSeparation(polymerArray, numUniqueMonomers, polymerLength, monomerID);
	var obj = chartData[0];
	keyList = Object.keys(obj);
	chart2.categoryField = keyList[0];
	chart2.graphs[0].valueField = keyList[1];
	chart2.dataProvider = chartData;
	chart2.categoryAxis.title = "Monomer " + monomerID + " Block Size";
	chart2.graphs[0].lineColor = colorArray[monomerID - 1];
	chart2.graphs[0].fillColors = colorArray[monomerID - 1];
	chart2.validateData();
	chart2.validateNow();
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