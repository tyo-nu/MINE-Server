

function mineDatabaseServices(url, auth, auth_cb) {

    this.url = url;
    var _url = url;
    var deprecationWarningSent = false;

    function deprecationWarning() {
        if (!deprecationWarningSent) {
            deprecationWarningSent = true;
            if (!window.console) return;
            console.log(
                "DEPRECATION WARNING: '*_async' method names will be removed",
                "in a future version. Please use the identical methods without",
                "the'_async' suffix.");
        }
    }

    var _auth = auth ? auth : { 'token' : '', 'user_id' : ''};
    var _auth_cb = auth_cb;


    this.quick_search = function (db, query, _callback, _errorCallback) {
    return json_call_ajax("mineDatabaseServices.quick_search",
        [db, query], 1, _callback, _errorCallback);
};

    this.quick_search_async = function (db, query, _callback, _error_callback) {
        deprecationWarning();
        return json_call_ajax("mineDatabaseServices.quick_search", [db, query], 1, _callback, _error_callback);
    };

    this.similarity_search = function (db, smiles, min_tc, _callback, _errorCallback) {
    return json_call_ajax("mineDatabaseServices.similarity_search",
        [db, smiles, min_tc], 1, _callback, _errorCallback);
};

    this.similarity_search_async = function (db, smiles, min_tc, _callback, _error_callback) {
        deprecationWarning();
        return json_call_ajax("mineDatabaseServices.similarity_search", [db, smiles, min_tc], 1, _callback, _error_callback);
    };

    this.database_query = function (params, _callback, _errorCallback) {
    return json_call_ajax("mineDatabaseServices.database_query",
        [params], 1, _callback, _errorCallback);
};

    this.database_query_async = function (params, _callback, _error_callback) {
        deprecationWarning();
        return json_call_ajax("mineDatabaseServices.database_query", [params], 1, _callback, _error_callback);
    };

    this.get_comps = function (db, ids, _callback, _errorCallback) {
    return json_call_ajax("mineDatabaseServices.get_comps",
        [db, ids], 1, _callback, _errorCallback);
};

    this.get_comps_async = function (db, ids, _callback, _error_callback) {
        deprecationWarning();
        return json_call_ajax("mineDatabaseServices.get_comps", [db, ids], 1, _callback, _error_callback);
    };

    this.get_rxns = function (db, ids, _callback, _errorCallback) {
    return json_call_ajax("mineDatabaseServices.get_rxns",
        [db, ids], 1, _callback, _errorCallback);
};

    this.get_rxns_async = function (db, ids, _callback, _error_callback) {
        deprecationWarning();
        return json_call_ajax("mineDatabaseServices.get_rxns", [db, ids], 1, _callback, _error_callback);
    };

    this.get_models = function (_callback, _errorCallback) {
    return json_call_ajax("mineDatabaseServices.get_models",
        [], 1, _callback, _errorCallback);
};

    this.get_models_async = function (_callback, _error_callback) {
        deprecationWarning();
        return json_call_ajax("mineDatabaseServices.get_models", [], 1, _callback, _error_callback);
    };

    this.get_adducts = function (_callback, _errorCallback) {
    return json_call_ajax("mineDatabaseServices.get_adducts",
        [], 1, _callback, _errorCallback);
};

    this.get_adducts_async = function (_callback, _error_callback) {
        deprecationWarning();
        return json_call_ajax("mineDatabaseServices.get_adducts", [], 1, _callback, _error_callback);
    };

    this.adduct_db_search = function (params, _callback, _errorCallback) {
    return json_call_ajax("mineDatabaseServices.adduct_db_search",
        [params], 1, _callback, _errorCallback);
};

    this.adduct_db_search_async = function (params, _callback, _error_callback) {
        deprecationWarning();
        return json_call_ajax("mineDatabaseServices.adduct_db_search", [params], 1, _callback, _error_callback);
    };

    this.pathway_search = function (pathway_query_params, _callback, _errorCallback) {
    return json_call_ajax("mineDatabaseServices.pathway_search",
        [pathway_query_params], 1, _callback, _errorCallback);
};

    this.pathway_search_async = function (pathway_query_params, _callback, _error_callback) {
        deprecationWarning();
        return json_call_ajax("mineDatabaseServices.pathway_search", [pathway_query_params], 1, _callback, _error_callback);
    };
 

    /*
     * JSON call using jQuery method.
     */
    function json_call_ajax(method, params, numRets, callback, errorCallback) {
        var deferred = $.Deferred();

        if (typeof callback === 'function') {
           deferred.done(callback);
        }

        if (typeof errorCallback === 'function') {
           deferred.fail(errorCallback);
        }

        var rpc = {
            params : params,
            method : method,
            version: "1.1",
            id: String(Math.random()).slice(2),
        };

        var beforeSend = null;
        var token = (_auth_cb && typeof _auth_cb === 'function') ? _auth_cb()
            : (_auth.token ? _auth.token : null);
        if (token != null) {
            beforeSend = function (xhr) {
                xhr.setRequestHeader("Authorization", token);
            }
        }

        var xhr = jQuery.ajax({
            url: _url,
            dataType: "text",
            type: 'POST',
            processData: false,
            data: JSON.stringify(rpc),
            beforeSend: beforeSend,
            success: function (data, status, xhr) {
                var result;
                try {
                    var resp = JSON.parse(data);
                    result = (numRets === 1 ? resp.result[0] : resp.result);
                } catch (err) {
                    deferred.reject({
                        status: 503,
                        error: err,
                        url: _url,
                        resp: data
                    });
                    return;
                }
                deferred.resolve(result);
            },
            error: function (xhr, textStatus, errorThrown) {
                var error;
                if (xhr.responseText) {
                    try {
                        var resp = JSON.parse(xhr.responseText);
                        error = resp.error;
                    } catch (err) { // Not JSON
                        error = "Unknown error - " + xhr.responseText;
                    }
                } else {
                    error = "Unknown Error";
                }
                deferred.reject({
                    status: 500,
                    error: error
                });
            }
        });

        var promise = deferred.promise();
        promise.xhr = xhr;
        return promise;
    }
}


