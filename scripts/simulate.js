numUniqueMonomers = 2;
currAnimation = 0;
///===============================================================================================================================
///'ENTER' KEY EVENT TRIGGER
///===============================================================================================================================
function checkEnter(event) {
	var keyPressed = event.keyCode;
	///console.log("keyPressed: ", keyPressed);
	if (keyPressed == 13) {
		simulate();
	}
}
//================================================================================================================================
//SIMULATION AND DYNAMIC VISUALIZATION
//================================================================================================================================
function simulate(demo = false) {
	colorArray = ["#FF6600", "#FCD202", "#B0DE09", "#0D8ECF", "#2A0CD0", 
	"#CD0D74", "#CC0000", "#00CC00", "#0000CC", "#DDDDDD", "#999999", "#333333", "#990000"];
	//Get all relevant inputs
	numUniqueMonomers = parseInt(document.getElementById("numUniqueMonomers").value);
	var totalNumMonomers = parseInt(document.getElementById("totalNumMonomers").value);
	var mRatio = parseInt(document.getElementById("mRatio").value);
	var conversion = parseInt(document.getElementById("conversion").value);
	var numRowsToShow = parseInt(document.getElementById("numRowsToShow").value);
	var graphTypeObj = document.getElementById("graph1Type");
	var graphType = graphTypeObj.options[graphTypeObj.selectedIndex].value;
	var animate = false;
	if (document.getElementById("animate").checked) {
		animate = true;
	}
	if (demo) {
		animate = true;
		numRowsToShow = 4;
		mRatio = 74;
		totalNumMonomers = 50000;
	}
	//Stops any current animation
	if (!demo) {
		stopAnimation();
	}

	//Initial variable calculations
	//console.log("Type of monomer1Ratio: ", typeof(monomer1Ratio));
	polymerLength = Math.floor(mRatio * conversion / 100);
	var numPolymers = Math.floor(totalNumMonomers / mRatio);
	var monomerRatioList = getMonomerRatios(numUniqueMonomers);
	console.log("monomerRatioList: ", monomerRatioList);

	//Set up visualizing canvas
	var canvas = document.getElementById("visual");
	var ctx = canvas.getContext("2d");
	var winWidth = window.innerWidth;
	var canvasWidth = parseInt(winWidth * 0.95);
	var squareLength = Math.min((canvasWidth - 10) / polymerLength, 15);
	canvas.setAttribute("width", canvasWidth);
	canvas.setAttribute("height", squareLength * numRowsToShow + 10);
	ctx.clearRect(0, 0, canvas.width, canvas.height);
	ctx.beginPath();
	ctx.translate(0.5, 0.5);
	ctx.lineWidth="1";
	ctx.strokeStyle = "black";
	//var rrList = [[monomer1RR, 1], [1, monomer2RR]];
	var monomerAmountsList = getMonomerAmounts(monomerRatioList, totalNumMonomers);
	initialMonomerAmountList = getMonomerAmounts(monomerRatioList, totalNumMonomers);
	//console.log("initialMonomerAmountList1: ", initialMonomerAmountList);
	rrList = getRateConstantRatios(numUniqueMonomers);
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
		} else {
		//If it is a 3-monomer or more system, use initial starting ratios to determing starting monomers
			var initChoices = [];
			var initWeightList = [];
			for (var monomerID = 1; monomerID <= numUniqueMonomers; monomerID++) {
				initChoices.push(monomerID);
				var initWeight = monomerAmountsList[monomerID - 1];
				initWeightList.push(initWeight);	
			}
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
	console.log("rrList 2: ", rrList);
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
	visualize(canvas, ctx, polymerArray, polymerLength, numRowsToShow, squareLength, animate, demo);
};

//================================================================================================================================
//FUNCTIONS FOR RETRIEVING USER INPUT DATA
//=================================================================================================================================

//Parses HTML document, returns list of monomer ratios
function getMonomerRatios(numUniqueMonomers) {
	var monomerRatioList = [];
	var inputTable = document.getElementById("inputTable");
	for (var index = 0; index < numUniqueMonomers; index++) {
		var row = inputTable.rows[index];
		var cell = row.cells[5];
		var inputObj = cell.children[0];
		var monomerRatio = parseFloat(inputObj.value); 
		monomerRatioList.push(monomerRatio);
	}
	return monomerRatioList;
};

//Parses HTML document, returns a list of lists of reactivity ratios

function getRateConstantRatios(numUniqueMonomers) {
	var rcList = [];
	var inputTable = document.getElementById("inputTable");
	if (numUniqueMonomers == 2) {
		var row1 = inputTable.rows[0];
		var cell1 = row1.cells[7];
		var inputObj1 = cell1.children[0];
		var rr1 = parseFloat(inputObj1.value);
		var singleRRList1 = [rr1, 1];
		var row2 = inputTable.rows[1];
		var cell2 = row2.cells[7];
		var inputObj2 = cell2.children[0];
		var rr2 = parseFloat(inputObj2.value);
		var singleRRList2 = [1, rr2];
		rcList = [singleRRList1, singleRRList2];
	} else {
		var rrList = [];
		for (var index = 0; index < numUniqueMonomers; index++) {
			var singleRRList = []; 
			rrList.push(singleRRList);
			for (var index2 = 0; index2 < numUniqueMonomers - 1; index2++) {
				var row = inputTable.rows[index];
				var cell = row.cells[7 + 2 * index2];
				var inputObj = cell.children[0];
				var reactivityRatio = parseFloat(inputObj.value);
				singleRRList.push(reactivityRatio);
			}
		}
		console.log("rrList: ", rrList);
		rcList = convertRRtoRC(numUniqueMonomers, rrList);
		//console.log("got here!");
	}
	console.log("rcList: ", rcList);
	return rcList;
};

//Input: a matrix of reactivity ratios
//Output: a larger matrix of relative rate constants derived from reactivity ratios

function convertRRtoRC(numUniqueMonomers, rrList) {
	var rcList = [];
	for (var rrIndex = 0; rrIndex < numUniqueMonomers; rrIndex++) {
		var singleRCList = [];
		var rrIndexToAccess = 0;
		var rrSubList = rrList[rrIndex];
		for (var rrIndex2 = 0; rrIndex2 < numUniqueMonomers; rrIndex2++) {
			if (rrIndex2 == rrIndex) {
				singleRCList.push(1);
			} else {
				rc =  1 / rrSubList[rrIndexToAccess];
				singleRCList.push(rc);
				rrIndexToAccess++;
			}
		}
		rcList.push(singleRCList); 
	}
	return rcList;
}

//==================================================================================================================================
//CREATING ITERATIVE INPUTS
//==================================================================================================================================

function createInputs(inputNum) {
	inputNum = parseFloat(inputNum);
	if (inputNum < 2) {
		return;
	}
	var inputTable = document.getElementById("inputTable");

	//Change HTML Objects for Monomer Ratios

	//Delete all current objects relating to monomer ratios

	for (var index = 0; index < numUniqueMonomers; index++) {
		var row = inputTable.rows[index];
		row.deleteCell(4);
		row.deleteCell(4);
	}

	//Iteratively create objects for Monomer Ratios

	for (var index = 0; index < inputNum; index++) {
		var labelElement = document.createElement("LABEL");;
		labelElement.innerHTML = "Monomer " + (index + 1) + " Ratio:";
		var inputElement = document.createElement("INPUT");
		inputElement.setAttribute("type", "number");
		inputElement.setAttribute("value", 1);
		var row = inputTable.rows[index];
		//console.log("row: ", "row"+index);
		var labelCell = row.insertCell(4);
		labelCell.appendChild(labelElement);
		var inputCell = row.insertCell(5);
		inputCell.appendChild(inputElement);
	}

	//Change HTML Objects for Monomer Reactivities

	//Delete all objects for monomer reactivities

	for (var index = 0; index < numUniqueMonomers; index++) {
			var row = inputTable.rows[index];
			if (numUniqueMonomers == 2) {
				row.deleteCell(-1);
				row.deleteCell(-1);
			} else {
				for (var index2 = 1; index2 <= numUniqueMonomers - 1; index2++) { 
					row.deleteCell(-1);
					row.deleteCell(-1);
				}
			}
	}

	//Iteratively create objects 	
	if (inputNum > 2) {
		for (var index = 0; index < inputNum; index++) {
			var row = inputTable.rows[index];
			for (var index2 = 0; index2 < inputNum; index2++) {
				if (index2 != index) {
					var labelElement = document.createElement("LABEL");
					var labelText = " rr" + (index + 1) + (index2 + 1);
					labelElement.innerHTML = labelText;
					var inputElement = document.createElement("INPUT");
					inputElement.setAttribute("type", "number");
					inputElement.setAttribute("value", 1);
					var labelCell = row.insertCell(-1);
					labelCell.appendChild(labelElement);
					var inputCell = row.insertCell(-1);
					inputCell.appendChild(inputElement);
				}
			}
		}
	}
	if (inputNum == 2) {
		for (var index = 0; index < 2; index++) {
			var row = inputTable.rows[index];
			var labelElement = document.createElement("LABEL");
			var labelText = "Reactivity Ratio" + (index + 1);
			labelElement.innerHTML = labelText;
			var inputElement = document.createElement("INPUT");
			inputElement.setAttribute("type", "number");
			inputElement.setAttribute("value", 1);
			var labelCell = row.insertCell(-1);
			labelCell.appendChild(labelElement);
			var inputCell = row.insertCell(-1);
			inputCell.appendChild(inputElement);
		}
	}

	//Update numUniqueMonomers

	numUniqueMonomers = inputNum;

};

//=========================================================================================================================================
//POLYMER ANALYSIS FUNCTIONS
//=========================================================================================================================================

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
};
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
};
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
};
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
};
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
};
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
};

//==========================================================================================================================================
//GRAPHING FUNCTIONS
//==========================================================================================================================================

function graphSetUp() {
	colorArray = ["#FF6600", "#FCD202", "#B0DE09", "#0D8ECF", "#2A0CD0", 
	"#CD0D74", "#CC0000", "#00CC00", "#0000CC", "#DDDDDD", "#999999", "#333333", "#990000"];
	// XY CHART
	chart = new AmCharts.AmXYChart();
	chart.startDuration = 0;
	chart["export"] = {"enabled": true};

	//X
	var xAxis = new AmCharts.ValueAxis();
	xAxis.title = "Monomer Index";
	xAxis.position = "bottom";
	xAxis.dashLength = 1;
	xAxis.axisAlpha = 0;
	xAxis.autoGridCount = true;
	xAxis.minimum = 0;
	xAxis.maximum = 100;
	chart.addValueAxis(xAxis);

	// Y
	var yAxis = new AmCharts.ValueAxis();
	yAxis.position = "left";
	yAxis.title = "Monomer Occurence";
	yAxis.dashLength = 1;
	yAxis.axisAlpha = 0;
	yAxis.autoGridCount = true;
	yAxis.minimum = 0;
	yAxis.maximum = 1;
	chart.addValueAxis(yAxis);

	//Remove all old graphs
	var graphList = chart.graphs;
	for (var i = 0; i < graphList.length; i++) {
		chart.removeGraph(graphList[i]);
	}
	//Add graphs depending on number of unique monomers
	for (var graphNumber = 1; graphNumber <= 2; graphNumber++) {
		var graph = new AmCharts.AmGraph();
	  	graph.lineColor = colorArray[graphNumber - 1];
		graph.balloonText = "x:[[x]] y:[[y]]";
		graph.xField = "x" + graphNumber;
		graph.yField = "y" + graphNumber;
		graph.lineAlpha = 0;
		graph.bullet = "triangleDown";
		graph.title = "Monomer" + graphNumber;
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
};

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
	chart2["export"] = {"enabled": true};

	//Value Axis
	var valueAxis = new AmCharts.ValueAxis();
	valueAxis.position = "bottom";
	valueAxis.minimum = 0;
	valueAxis.title = "Normalized Counts";
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
};

function formatXYGraphs() {
	var chartData = getMonomerComposition(polymerArray, numUniqueMonomers, polymerLength);

	//Remove all old graphs
	var graphList = chart.graphs;
	console.log("graphList.length", graphList.length);
	initGraphListLen = graphList.length;
	for (var i = 0; i < initGraphListLen; i++) {
		console.log("graphList: ", graphList, i);
		chart.removeGraph(graphList[0]);
	}
	console.log("lengthafter: ", graphList.length);
	console.log("numUniqueMonomers: ", numUniqueMonomers);
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
	//chart.removeGraph(graphList[0]);
};
function formatHistGraph(monomerID) {
	if (monomerID > numUniqueMonomers) {
		return;
	}
	monomerID = parseInt(monomerID);
	console.log("monomerID: ", monomerID);
	var chartData = getMonomerSeparation(polymerArray, numUniqueMonomers, polymerLength, monomerID);
	var obj = chartData[0];
	keyList = Object.keys(obj);
	chart2.categoryField = keyList[0];
	chart2.graphs[0].valueField = keyList[1];
	chart2.dataProvider = chartData;
	chart2.categoryAxis.title = "Monomer " + monomerID + " Run Length";
	chart2.graphs[0].lineColor = colorArray[monomerID - 1];
	chart2.graphs[0].fillColors = colorArray[monomerID - 1];
	chart2.validateData();
	chart2.validateNow();
};
function setGraph(type) {
	//console.log("type: ", type);
	formatXYGraphs();
	var chartData;
	switch (type) {
		case "Monomer Occurences":
			chartData = getMonomerComposition(polymerArray, numUniqueMonomers, polymerLength);
			//console.log("chartData: ", chartData);
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
		case "Run length":
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
	//console.log("got here");
	chart.invalidateSize();
};

//=================================================================================================================================
//VISUALIZING FUNCTIONS
//=================================================================================================================================
function visualize(canvas, ctx, polymerArray, polymerLength, numRowsToShow, length, animate, demo) {
	if (currAnimation) {
		cancelAnimationFrame(currAnimation);
	}
	ctx.clearRect(0, 0, canvas.width, canvas.height);
	ctx.beginPath();
	var blockList = [];
	for (var len = 0; len < polymerLength; len++) {
		for (var polymerIndex = 0; polymerIndex < numRowsToShow; polymerIndex ++) {
			monomerID = polymerArray[polymerIndex][len];
			var color = colorArray[monomerID - 1];
			var block = {
				bcolor: color,
				ulX: len * length,
				ulY: polymerIndex * length
			}
			blockList.push(block);
			if (!animate) {
				draw();
				function draw() {
					ctx.beginPath()
					ctx.fillStyle = color;
					ctx.rect(len*length, polymerIndex*length, length, length);
					ctx.fill();
					ctx.stroke();
				}
			}

		}
	}
	var lastIndex = 0;
	function draw() {
		var fadeDelay =  10;
		ctx.clearRect(0, 0, canvas.width, canvas.height);
		for (var index = 0; index < lastIndex; index++) {
			ctx.globalAlpha = 1;
			var initFade = 0.2;
			var maxFade = 1;
			if (index > lastIndex - fadeDelay * numRowsToShow) {
				fadeStep = (maxFade - initFade) / (fadeDelay * numRowsToShow);
				fadeDistance = Math.abs(lastIndex - index);
				fade = initFade + fadeStep * fadeDistance;
				ctx.globalAlpha = fade;
			}
			if (index < blockList.length){ 
				block = blockList[index];
				ctx.beginPath()
				ctx.fillStyle = block.bcolor;
				ctx.rect(block.ulX, block.ulY, length, length);
				ctx.fill();
				ctx.stroke();
			}
		}
		if (lastIndex > blockList.length + fadeDelay * numRowsToShow) {
			//animation.cancel();
			// for (var index = 0; index < lastIndex; index++) {
			// 	ctx.globalAlpha = 1;
			// 	block = blockList[index];
			// 	ctx.beginPath()
			// 	ctx.fillStyle = block.bcolor;
			// 	ctx.rect(block.ulX, block.ulY, length, length);
			// 	ctx.fill();
			// 	ctx.stroke();
			// }
		}
		if (lastIndex < blockList.length + fadeDelay * numRowsToShow) {
			lastIndex += 1;
		}
		currAnimation = requestAnimationFrame(draw)
	}
	if (animate) {
		draw();
	}
	// var ctx = canvas.getContext("2d");
	// ctx.fillRect(0,0,50,50);
	// ctx.fillRect(100,0,50,50);
	// ctx.fillRect(200,0,50,50);
}
function stopAnimation() {
	interrupt = true;
}
function setCanvasScalingFactor() {
	return window.devicePixelRatio || 1;
}
//================================================================================================================================
//UTILITY FUNCTIONS
//================================================================================================================================
//Creates an array from [1, 2, 3,......, N]
function createRangeArray(length) {
	var array = [];
	for (var i = 1; i <= length; i++) {
		array.push(i);
	}
	return array;
};
//Sums the numerical contents of an array
function sumArray(array) {
	return array.reduce(sum, 0)
};
function sum(a, b) {
	return a + b;
};
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
//sleep function
function sleep(milliseconds) {
  var start = new Date().getTime();
  for (var i = 0; i < 1e7; i++) {
    if ((new Date().getTime() - start) > milliseconds){
      break;
    }
  }
}
//========================================================================================================================================
//DATA
//========================================================================================================================================
