<?php
// open config-file
require("php_config.php");

$connection = new mysqli("$host", "$username", "$password", $db_name);
if ($connection->connect_errno)
    die("cannot connect");

// get values from database
$res = $connection->query("SELECT aminoacid FROM substrates");
$aminoacids = array();

while ($row = $res->fetch_row()) {
    array_push($aminoacids, $row[0]);
}

//var_dump($aminoacids);


// insert modifications of aminoacids in json
//while ($row = $res->fetch_row())
//$arr = array ('a'=>1,'b'=>2,'c'=>3,'d'=>4,'e'=>5);

//  echo json_encode($arr);

// test with array
//$aminoacids = array("Alanin", "Ornithin", "Histidin");

?>

<html>
    <head>
        <title>NRPS Designer</title>
        <style type="text/css">
              body
        {
            
            background-image: radial-gradient(ellipse farthest-corner at left top, #FFFFFF 0%, #000A0F 100%);
            /*background-repeat: no-repeat;*/
            height: 100%;
            
        }
        </style>
        <script type="text/javascript">
            function makeDropdown(obj) {
                var p = document.getElementById("aminoacids");
                if (obj.value > 100) {
                    obj.value = 100;
                    alert("Maximal 100 amino acids allowed");
                }
                var inner = "";
                for (var i = 0; i < obj.value; ++i) {
<?php
echo "inner += \"<select name='aminoacids[]' size='1'>";
foreach ($aminoacids as $aminoacid) {
    echo "<option>$aminoacid</option>";
    
}
echo "</select><input type='radio' name='aa[\" + i + \"]' value='L'> L-conf.<br><input type='radio' name='aa[]' value='D'> D-conf.<br>\"";
?>
                    ;p.innerHTML = inner;
                }
            }
        </script>
    </head>
    <body>
<?php
if (!isset($_POST['aminoacids']) || !isset($_POST['aa']) || count($_POST['aminoacids']) != count($_POST['aa'])) {
?>
        <form action="<?php echo $_SERVER['PHP_SELF'];?>" method="post">
            organism: <input name="organism" type="text" size"30" maxlength="30";>
            length: <input name="length" type="text" size="30" maxlength="30" onkeyup="makeDropdown(this);">
            <p id="aminoacids">
                
            </p>
            <input type="submit" value="submit">
        </form>
<?php
} else {
    $oXMLWriter = new XMLWriter;
    $oXMLWriter->openMemory();
    $oXMLWriter->setIndent(true);
    $oXMLWriter->setIndentString("    ");
    $oXMLWriter->startDocument('1.0', 'UTF-8');
    
    $oXMLWriter->startElement('nrp');
        
            for ($i = 0; $i < count($_POST['aminoacids']); ++$i) {
                $oXMLWriter->startElement('monomer');
                $oXMLWriter->startElement('name');
                $oXMLWriter->text($_POST['aminoacids'][$i]);
                $oXMLWriter->endElement();
                
                $oXMLWriter->startElement('configuration');
                $oXMLWriter->text($_POST['aa'][$i]);
                $oXMLWriter->endElement();
                
                //$oXMLWriter->startElement('modification');
                //$oXMLWriter->endElement();
                
        $oXMLWriter->endElement();
            }
    $oXMLWriter->endElement();
    $oXMLWriter->endDocument();
    echo $oXMLWriter->outputMemory(TRUE);
}
?>
    </body>
</html>
