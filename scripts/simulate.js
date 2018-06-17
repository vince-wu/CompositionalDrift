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
	var rrList = [monomer1RR, monomer2RR];
	var monomerAmountsList = getMonomerAmounts(monomerRatioList, totalNumMonomers);
	console.log(monomerAmountsList);
	var polymerArray = initiation(numUniqueMonomers, numPolymers, rrList, monomerAmountsList);
	console.log(polymerArray);

}
function initiation(numUniqueMonomers, numPolymers, rrList, monomerAmountsList) {
	var polymerArray = [];
	var currNumPolymers;
	for (currNumPolymers = 0; currNumPolymers < numPolymers; currNumPolymers++) {
		var spec = {};
		//if it is a 2-monomer system, use the Mayo Lewis Equation to Determine Starting Monomers
		//Weights that determine probabilty of monomer starting the chain are calculated
		if (numUniqueMonomers == 2) {
			var f1 = monomerAmountsList[0];
			var f2 = monomerAmountsList[1];
			var r1 = rrList[0];
			var r2 = rrList[1];
			var weight = (r1*f1**2 + f1*f2) / (r1*f1**2 + 2*f1*f2 + r2*f2**2);
			spec[1] = weight;
			spec[2] = 1 - weight;
		} 
		//Use a weighted random selector to choose the starting monomer
		var startingMonomer = weightedRand(spec);
		console.log("startingMonomer: ", startingMonomer);
		//Create a new polymer chain with the selected starting monomer
		polymerArray.push([startingMonomer]);
		//Remove one of that monomer from the pool
		monomerAmountsList[startingMonomer - 1] --;
	}
	return polymerArray;

}
//Gives the starting number of each unique monomer in list form
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
//Sums the numerical contents of an array
function sumArray(array) {
	return array.reduce(sum, 0)
}
function sum(a, b) {
	return a + b;
}
function weightedRand(spec) {
  var i, sum=0, r=Math.random();
  for (i in spec) {
    sum += spec[i];
    if (r <= sum) return i;
  }
}