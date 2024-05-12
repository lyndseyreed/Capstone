# Spring 2024 Bioinformatics Capstone Project
## Created by Lyndsey Reed


## **General Overview**
This project looks to take in protein sequence information and display a 3D model. Due to the novel nature of the input proteins, the structure must be predicted. This code will utilize ESMFold as the main predictive resource for the models and Streamlit web based application as the medium for the output visual. This README provides guidelines on how to install, configure, and use this tool.

## **Features**

* **ESMFold Integration:** Predicts protein structures from fasta files.
* **Streamlit Visualization:** Interactive 3D visualization of the predicted protein structures.
* **File Input/Output:** Accepts fasta files as input and outputs PDB files and 3D structure.

## **Installation**
### **Requirements:**
* Python 3.8 or later
* pip for package management
* 20-50GB free RAM for computations

### **Steps:**

1. **Clone the Repository:**
``` python
git clone https://github.com/lyndseyreed/Capstone.git
cd Capstone
```
2. **Load the Requirements:**
```python
pip install -r requirements.txt
```

## **Usage**
1. Uplaod your FASTA format files into the **Input** folder.

2. To extract the sequences run:
```python
python3 sequence_extractor.py
```
3. To get the model predictions run:
```pyython3
python3 prediction_generator_binary.py
```
4. To prepare the file for 3D rendering, you must remove the binary by running:
```python
python3 remove_binary_to_pdb.py
```
5. To view the rendered files in the web based application run:
```python
streamlit run model_generator.py
```

**11MAY24: Please note that this project is a work in progress, there are additional files that are not currently in use but will be deployed at a later date. The model prediction currently takes an extended period of time and computational power and will be improved or replaced.**

## Third Party Code Citation
Function **convert_to_pdb** in file **prediction_generator_binary.py** was sourced from https://github.com/facebookresearch/esm.git under the MIT License:

>MIT License
>
>Copyright (c) Meta Platforms, >Inc. and affiliates.
>
>Permission is hereby granted, >free of charge, to any person >obtaining a copy
>of this software and associated >documentation files (the >"Software"), to deal
>in the Software without >restriction, including without limitation the rights
>to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
>copies of the Software, and to permit persons to whom the Software is
>furnished to do so, subject to the following conditions:
>
>The above copyright notice and this permission notice shall be included in all
>copies or substantial portions of the Software.
>
>THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
>IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
>FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
>AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
>LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
>OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
>SOFTWARE.


