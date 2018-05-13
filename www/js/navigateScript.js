$(document).ready(function () {

    loadPage();

    $('a[href="#"]').click(function () {
        window.history.pushState('', null, '/#home');
        loadHome();
    });

    $('a[href="#home"]').click(function () {
        window.history.pushState('', null, '/#home');
        loadHome();
    });

    $('a[href="#project"]').click(function () {
        window.history.pushState('', null, '/#project');
        loadProject();
    });

    function loadPage() {
        switch (window.location.hash) {
            case "#project":
                loadProject();
            default:
                loadHome();
        }
    }

    function loadProject() {
        $('#page-content').load("project.html #page-content > *", function (responseText, textStatus, XMLHttpRequest) {
            loadNetwork();
        });
    }

    function loadNetwork(){
        $('#network').load("network.html > *", function (responseText, textStatus, XMLHttpRequest) {
            jQuery.getScript("js/networkScript.js");
        });
    }

    function loadHome() {
        window.history.pushState('', null, '/#home');
        $('#page-content').load("home.html #page-content > *", function (responseText, textStatus, XMLHttpRequest) {
            jQuery.getScript("js/homeScript.js");
            $('#workflow').load("workflow/index.html > *", function (responseText, textStatus, XMLHttpRequest) {
            });
        });
    }
});

