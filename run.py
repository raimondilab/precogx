from flask import Flask, render_template, request
import json

app = Flask(__name__)

@app.route('/output', methods=['GET', 'POST'])
def output():
    #data = request.form['input']
    if request.method == 'POST':
        #print (request.form['input'])
        path = "static/OL820/out.json"
        return render_template('embedded.html', name='name', path_to_json_output=json.dumps(path))
    return ("<h3>It was a GET request</h3>")

@app.route('/', methods=['GET', 'POST'])
@app.route('/home', methods=['GET', 'POST'])
def home():
    return render_template('home.html')

if __name__ == '__main__':
    app.run(debug=True, host='0.0.0.0')
