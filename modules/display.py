from PyQt5 import QtWidgets, QtCore
from pkg_resources import parse_version
import requests
from requests.exceptions import HTTPError
import json
import re

def initDisplay(self):
	#Set initial display background color
	self.displayView.setStyleSheet("background-color: black; padding: 50px;")
	self.displayView.setAlignment(QtCore.Qt.AlignCenter)
	self.displayView.setTextInteractionFlags(QtCore.Qt.LinksAccessibleByMouse)

	versionText = getVersionText(self)

	#self.displayView.setBackgroundBrush(brush)
	txt = '\
	<div style="color: #FFF;">\
		<h2 style="text-align: center">\
			Welcome to Compositional Drift!\
		</h2>\
		<span style="text-align: center;font-size:9pt;">\
			<br>\
			Author: Vincent Wu (vincent.wu@berkeley.edu)\
			<br>\
			Contributors: Anton A. A. Smith, Aaron Hall\
		</span>\
		<p style="text-align:left;font-size:9pt;">\
			This program uses the Mayo-Lewis Equation and Monte Carlo Methods to simulate controlled polymer\
			reactions. Fill out your desired reaction parameters and click "Simulate" or press the "enter" key.\
			For more detailed instructions and to access the open-source code,\
	 		please refer to the project'"'"'s github page: \
	 		<a style="color: #FF0;" href="https://github.com/vince-wu/CompositionalDrift">\
			https://github.com/vince-wu/CompositionalDrift </a>.\
		</p>\
		 	{}\
	</div>'.format(versionText)
	text = self.displayView.setText(txt)
	#text.setDefaultTextColor(QtCore.Qt.yellow)
	self.displayView.setOpenExternalLinks(True)
	#text.setReadOnly(False)

def rr_initDisplay(self):
	#Set initial display background color
	self.rr_textBrowser.setStyleSheet("background-color: black; padding: 50px;")
	self.rr_textBrowser.setAlignment(QtCore.Qt.AlignCenter)
	self.rr_textBrowser.setTextInteractionFlags(QtCore.Qt.LinksAccessibleByMouse)

	versionText = getVersionText(self)

	#self.displayView.setBackgroundBrush(brush)
	txt = '\
	<div style="color: #FFF;">\
		<h2 style="text-align: center">\
			Welcome to Reactivity Ratio Calculator!\
		</h2>\
		<span style="text-align: center;font-size:9pt;">\
		</span>\
		<p style="text-align:left;font-size:9pt;">\
			This program estimates reactivity ratios for 2-monomer systems from experimental data.\
			For more detailed instructions and to access the open-source code,\
	 		please refer to the project'"'"'s github page: \
	 		<a style="color: #FF0;" href="https://github.com/vince-wu/CompositionalDrift">\
			https://github.com/vince-wu/CompositionalDrift </a>.\
		</p>\
		 	{}\
	</div>'.format(versionText)

	text = self.rr_textBrowser.setText(txt)
	#text.setDefaultTextColor(QtCore.Qt.yellow)
	self.rr_textBrowser.setOpenExternalLinks(True)
	#text.setReadOnly(False)


def getVersionText(self):
	versionText = ''
	debug = False
	connected = True
	#self.version = "v1.8"
	if debug:
		newestVersion = 'v2.0.0'
	try:
		bodyText = ''
		if not debug:
			#Get program's newest version from Github API, compare to current version
			json_response = getJson()
			if json_response:
				newestVersion = json_response['tag_name']
				bodyText = json_response['body']
				bodyText = re.sub(r'(<li>)', '<li>- ', bodyText)
				bodyText = re.sub(r'(Changelog:)', '{} Changelog:'.format(newestVersion), bodyText)
				newestVersionNumber = re.sub('[^0-9,.]','',newestVersion)
			else:
				connected = False


		currentVersionNumber = re.sub('[^0-9,.]','',self.version)
		#Alert user based on if version is updated, ahead, or behind
		if connected == False:
			versionText = '\
			<p style="color: #FFF; font-size:9pt; text-align:center">\
				<br>\
				<span style="color: #FF0">\
					***\
				</span>\
				To see if this program is up to date, please connect to the internet or \
				visit \
				<a style="color: #FF0" href="https://github.com/vince-wu/CompositionalDrift/releases/latest">\
					this link\
				</a>\
				.\
				<span style="color: #FF0">\
					***\
				</span>\
			</p>\
			'
		elif parse_version(newestVersionNumber) > parse_version(currentVersionNumber):
			versionText = '\
			<p style="color: #FFF; font-size:9pt; text-align:center">\
				<br>\
				<span style="color: #FF0">\
					***\
				</span>\
				An updated version of this software ({}) \
				is available\
				<a style="color: #FF0" href="https://github.com/vince-wu/CompositionalDrift/releases/latest">\
					here\
				</a>\
				!\
				<span style="color: #FF0">\
					***\
				</span>\
			</p>\
			<p style="color: #FFF; font-size:8pt">\
				{}\
			</p>\
			'.format(newestVersion, bodyText)

		elif parse_version(newestVersionNumber) == parse_version(currentVersionNumber):
			versionText = '\
			<p style="color: #FFF; font-size:9pt; text-align:center">\
				<br>\
				<span style="color: #FF0">\
					***\
				</span>\
					This program version is up to date!\
				<span style="color: #FF0">\
					***\
				</span>\
				<p style="color: #FFF; font-size:8pt">\
					{}\
				</p>\
			</p>\
			'.format(bodyText)
		elif parse_version(newestVersionNumber) < parse_version(currentVersionNumber):
			versionText = '\
			<p style="color: #FFF; font-size:9pt; text-align:center">\
				<br>\
				<span style="color: #FF0">\
					***\
				</span>\
					You are running an unpublished pre-release version of the program!  \
					Please report any bugs that you find to vincent.wu@berkeley.edu.\
				<span style="color: #FF0">\
					***\
				</span>\
				<p style="color: #FFF; font-size:8pt">\
					{}\
				</p>\
			</p>\
			'.format(bodyText)

	except Exception as e:
		pass

	return versionText

def getJson():
	try:
		response = requests.get(
			'https://api.github.com/repos/vince-wu/CompositionalDrift/releases/latest'
			)
		json_response = response.json()
		return json_response
	except HTTPError as http_err:
		return None
	except Exception as err:
		return None
